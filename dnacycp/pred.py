import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import tensorflow as tf
import os
from typing import List, Tuple

network_final_original = keras.models.load_model("irlstm")
detrend_int_original = 0.029905550181865692
detrend_slope_original = 0.973293125629425
normal_mean_original = -0.18574825868055558
normal_std_original = 0.4879013326394626

network_final_smooth = keras.models.load_model("irlstm_smooth")
detrend_int_smooth = 0.001641373848542571
detrend_slope_smooth = 1.0158132314682007
# Mean and stdev of smoothed C0 for Tiling library:
# (calculated from/on the scale of normalized Cn values)
normal_mean_smooth = -0.011196041799376931
normal_std_smooth = 0.651684644408004

def dnaOneHot(sequence):
    code = np.array([
        [1, 0, 0, 0],  # A / a
        [0, 1, 0, 0],  # C / c
        [0, 0, 1, 0],  # G / g
        [0, 0, 0, 1],  # T / t
        [0, 0, 0, 0]   # N / n
    ])
    
    mapping = np.zeros(128, dtype=int)
    mapping[ord('A')] = 0
    mapping[ord('C')] = 1
    mapping[ord('G')] = 2
    mapping[ord('T')] = 3
    mapping[ord('N')] = 4
    mapping[ord('a')] = 0
    mapping[ord('c')] = 1
    mapping[ord('g')] = 2
    mapping[ord('t')] = 3
    mapping[ord('n')] = 4

    indices = np.fromiter((mapping[ord(char)] for char in sequence), dtype=int)
    onehot_encoded_seq = code[indices]

    return onehot_encoded_seq

def construct_predict_fn(model):
    @tf.function(reduce_retracing=True)
    def predict_fn(x):
        return model(x, training=False)
    return predict_fn

def cycle_fasta(inputfile:str, outputbase:str, smooth:bool=True, chunk_size=None, num_threads=None):
    """
    Make predictions for a given FASTA file.

    Parameters
    ----------
    inputfile : str
        The path to the FASTA file to predict.
    outputbase : str
        The base name of the output files.
    smooth : bool, optional
        Whether to use the smoothed model or not. The default is True.
        smooth=True corresponds to DNAcycP2, smooth=False corresponds to DNAcycP

    Notes
    -----
    The output files will be named as `<outputbase>_cycle_<chrom>.txt`, where `<chrom>` is the chromosome ID from the FASTA file.
    """
    genome_file = SeqIO.parse(open(inputfile),'fasta')

    network_final = network_final_smooth if smooth else network_final_original
    detrend_int = detrend_int_smooth if smooth else detrend_int_original
    detrend_slope = detrend_slope_smooth if smooth else detrend_slope_original
    normal_mean = normal_mean_smooth if smooth else normal_mean_original
    normal_std = normal_std_smooth if smooth else normal_std_original

    if chunk_size is None:
        chunk_size = 100000
    if num_threads is None:
        num_threads = os.cpu_count()

    predict_fn = construct_predict_fn(network_final)

    if smooth:
        print(f"Making smooth predictions (DNAcycP2)\n\n", flush=True)
    else:
        print(f"Making predictions (DNAcycP)\n\n", flush=True)

    def process_chunk(args) -> Tuple[np.ndarray, np.ndarray]:
        """
        Process a chunk of the sequence data with both forward and reverse predictions.
        """
        onehot_sequence, start_ind, end_ind, network_final = args
        
        # Extract local sequence window
        ind_local = np.arange(start_ind, end_ind)
        onehot_sequence_local = onehot_sequence[np.arange(ind_local[0] - 25, ind_local[-1] - 24)[:, None] + np.arange(50)]
        onehot_sequence_local = onehot_sequence_local.reshape((-1, 50, 4, 1))
        
        # Create reverse sequence
        onehot_sequence_local_reverse = np.flip(onehot_sequence_local, [1, 2])
        
        # Convert to TensorFlow tensors with fixed shape
        onehot_sequence_local = tf.cast(onehot_sequence_local, tf.float32)
        onehot_sequence_local_reverse = tf.cast(onehot_sequence_local_reverse, tf.float32)
        
        # Make predictions using the optimized prediction function
        fit_local = predict_fn(onehot_sequence_local).numpy().reshape(-1)
        fit_local_reverse = predict_fn(onehot_sequence_local_reverse).numpy().reshape(-1)
        
        return fit_local, fit_local_reverse
    
    for fasta in genome_file:
        chrom = fasta.id
        genome_sequence = str(fasta.seq)
        print(f"Sequence length for ID {chrom}: {len(genome_sequence)}", flush=True)
        onehot_sequence = dnaOneHot(genome_sequence)
        onehot_sequence = array(onehot_sequence)
        onehot_sequence = onehot_sequence.reshape((onehot_sequence.shape[0],4,1))
        print("Predicting cyclizability...", flush=True)

        print(f"Chunk size: {chunk_size}, num threads: {num_threads}", flush=True)

        sequence_length = onehot_sequence.shape[0] - 49
        start_indices = range(25, sequence_length + 25, chunk_size)
        end_indices = [min(start + chunk_size, sequence_length + 25) for start in start_indices]
        
        chunk_args = [
            (onehot_sequence, start, end, predict_fn) 
            for start, end in zip(start_indices, end_indices)
        ]
        
        fit = []
        fit_reverse = []
        
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            for i, (fit_local, fit_local_reverse) in enumerate(executor.map(process_chunk, chunk_args)):
                fit.append(fit_local)
                fit_reverse.append(fit_local_reverse)
                
                # Progress reporting for long sequences
                sequences_processed = len(fit_local)
                total_processed = sum(len(chunk) for chunk in fit)
                if onehot_sequence.shape[0] > 10**7:
                    print(f"\t Completed predictions on {total_processed} out of {sequence_length} sequences", flush=True)

        fit = np.concatenate(fit)  # Assuming fit is a list of arrays
        fit_reverse = np.concatenate(fit_reverse)
        fit = detrend_int + (fit + fit_reverse) * detrend_slope / 2
        fit2 = fit * normal_std + normal_mean
        n = fit.shape[0]
        positions = np.arange(25, 25 + n)
        fitall = np.column_stack((positions, fit, fit2))
        fitall = pd.DataFrame(fitall, columns=["position", "c_score_norm", "c_score_unnorm"])
        fitall = fitall.astype({"position": int})
        fitall.to_csv(outputbase+"_cycle_"+chrom+".txt", index = False)
        print("Output file: "+outputbase+"_cycle_"+chrom+".txt", flush=True)

def cycle_txt(inputfile:str, outputbase:str, smooth:bool=True):
    """
    Make predictions for a given TXT file.

    Parameters
    ----------
    inputfile : str
        The path to the TXT file to predict.
    outputbase : str
        The base name of the output files.
    smooth : bool, optional
        Whether to use the smoothed model or not. The default is True.
        smooth=True corresponds to DNAcycP2, smooth=False corresponds to DNAcycP

    Notes
    -----
    The output files will be named as `<outputbase>_cycle_norm.txt` and `<outputbase>_cycle_unnorm.txt`, where `<outputbase>` is the base name given as an argument.
    """
    network_final = network_final_smooth if smooth else network_final_original
    detrend_int = detrend_int_smooth if smooth else detrend_int_original
    detrend_slope = detrend_slope_smooth if smooth else detrend_slope_original
    normal_mean = normal_mean_smooth if smooth else normal_mean_original
    normal_std = normal_std_smooth if smooth else normal_std_original

    if smooth:
        print(f"Making smooth predictions (DNAcycP2)\n\n")
    else:
        print(f"Making predictions (DNAcycP)\n\n")

    with open(inputfile) as f:
            input_sequence = f.readlines()
    X = []
    all50 = True
    print("Reading sequences...")
    for loop_sequence in input_sequence:
        loop_sequence = loop_sequence.rstrip()
        if len(loop_sequence) != 50:
            all50=False
        X.append(dnaOneHot(loop_sequence))
    if all50:
        print("Predicting cyclizability...")
        X = array(X)
        X = X.reshape((X.shape[0],50,4,1))
        X_reverse = np.flip(X,[1,2])

        model_pred = network_final.predict(X)
        model_pred_reverse = network_final.predict(X_reverse)

        model_pred = detrend_int + (model_pred + model_pred_reverse) * detrend_slope/2
        output_cycle = model_pred.flatten()
        output_cycle2 = np.array([item * normal_std + normal_mean for item in output_cycle])
        output_cycle = list(output_cycle)
        output_cycle2 = list(output_cycle2)
    else:
        print("Not all sequences are length 50, predicting every subsequence...")
        output_cycle = []
        lenX = len(X)
        for j, onehot_loop in enumerate(X):
            l = len(onehot_loop)
            onehot_loop = array(onehot_loop)
            onehot_loop = onehot_loop.reshape((l,4,1))
            onehot_loops = []
            for i in range(l-49):
                onehot_loops.append(onehot_loop[i:i+50])
            onehot_loops = array(onehot_loops)
            onehot_loops_reverse = np.flip(onehot_loops,[1,2])
            if l > 1000:
                # Provide status bar for long sequences:
                cycle_local = network_final.predict(onehot_loops)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse)
            else:
                # No status bar for short sequences (verbose=0):
                cycle_local = network_final.predict(onehot_loops, verbose=0)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse, verbose=0)
            cycle_local = detrend_int + (cycle_local + cycle_local_reverse) * detrend_slope/2
            cycle_local = cycle_local.reshape(cycle_local.shape[0])
            output_cycle.append(cycle_local)
            if j%10==9:
                print(f"Completed {j+1} out of {lenX} total sequences")
        output_cycle2 = [item * normal_std + normal_mean for item in output_cycle]
    with open(outputbase+"_cycle_norm.txt", "w") as file:
        for row in output_cycle:
            if isinstance(row, (np.floating, float)):
                s = str(row)
            else:
                s = " ".join(map(str, row))
            file.write(s+'\n')
    with open(outputbase+"_cycle_unnorm.txt", "w") as file:
        for row in output_cycle2:
            if isinstance(row, (np.floating, float)):
                s = str(row)
            else:
                s = " ".join(map(str, row))
            file.write(s+'\n')
