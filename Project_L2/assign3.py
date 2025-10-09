# Design an application using AI, which contains a GUI that allows the user to select a FASTA file.
# The content of the file should be analyzed by using a sliding window of 30 positions.
# The content of each sliding window should be used in order to extract the relative frequencies of the symbols found in the alphabet of the sequence
# Thus, your input should be the DNA sequence from the FASTA file, and the output should be values of the relative frequencies for each symbol from the alphabet (sequence).
# Translate in lines on a chart, thus your chart in the case of DNA should contain 4 lines, one for each symbol found over the sequence.

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio import SeqIO
import numpy as np

DNA_ALPHABET = ['A', 'C', 'G', 'T']
WINDOW_SIZE = 30

class FastaAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Nucleotide Frequency Analyzer")
        self.root.geometry("1000x750") 

        self.current_file_path = None 

        top_frame = tk.Frame(self.root, pady=5)
        top_frame.pack(side="top", fill="x")
        
        control_frame = tk.Frame(self.root, pady=5)
        control_frame.pack(side="top", fill="x")

        chart_frame = tk.Frame(self.root)
        chart_frame.pack(side="top", fill="both", expand=True, padx=10, pady=10)

        self.btn_select_file = tk.Button(
            top_frame,
            text="Select FASTA File",
            command=self.select_file,
            font=("Helvetica", 12)
        )
        self.btn_select_file.pack(side="left", padx=20)

        self.lbl_file_path = tk.Label(
            top_frame,
            text="No file selected",
            font=("Helvetica", 10),
            fg="gray"
        )
        self.lbl_file_path.pack(side="left")

        self.smoothing_enabled = tk.BooleanVar(value=True)
        self.smoothing_window_size = tk.IntVar(value=51) # Default smoothing window

        lbl_smoothing = tk.Label(control_frame, text="Smoothing Controls:", font=("Helvetica", 10, "bold"))
        lbl_smoothing.pack(side="left", padx=(20, 10))

        self.chk_enable_smoothing = tk.Checkbutton(
            control_frame,
            text="Enable Smoothing",
            variable=self.smoothing_enabled,
            command=self._trigger_replot
        )
        self.chk_enable_smoothing.pack(side="left")

        lbl_window = tk.Label(control_frame, text="Strength (Window Size):")
        lbl_window.pack(side="left", padx=(20, 5))

        self.spin_smoothing_window = tk.Spinbox(
            control_frame,
            from_=3,
            to=201,
            increment=2, 
            textvariable=self.smoothing_window_size,
            width=5,
            command=self._trigger_replot
        )
        self.spin_smoothing_window.pack(side="left")

        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        
        self.setup_initial_chart()

    def setup_initial_chart(self):
        self.ax.clear()
        self.ax.set_title("Sliding Window Nucleotide Frequencies")
        self.ax.set_xlabel("Window Start Position")
        self.ax.set_ylabel("Relative Frequency")
        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.ax.text(0.5, 0.5, "Please select a FASTA file to begin analysis.",
                     horizontalalignment='center', verticalalignment='center',
                     transform=self.ax.transAxes, fontsize=14, color='gray')
        self.canvas.draw()

    def select_file(self):
        file_path = filedialog.askopenfilename(
            title="Select a FASTA file",
            filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*"))
        )
        if file_path:
            self.current_file_path = file_path # Store the file path
            self.lbl_file_path.config(text=file_path, fg="black")
            self.process_and_plot()

    def _trigger_replot(self):
        if self.current_file_path:
            self.process_and_plot()

    def apply_smoothing(self, data, window_size):
       
        if window_size <= 1:
            return np.array(data)
        
        kernel = np.ones(window_size) / window_size
        smoothed_data = np.convolve(data, kernel, mode='valid')
        return smoothed_data

    def process_and_plot(self):

        if not self.current_file_path:
            return

        try:
            fasta_sequences = SeqIO.parse(open(self.current_file_path), 'fasta')
            first_record = next(fasta_sequences)
            sequence = str(first_record.seq).upper()
            
            if len(sequence) < WINDOW_SIZE:
                messagebox.showerror("Error", f"Sequence length is too short.")
                return

            frequencies, window_positions = self.sliding_window_analysis(sequence)

            is_smoothed = self.smoothing_enabled.get()
            smooth_win_size = self.smoothing_window_size.get()

            if is_smoothed and len(window_positions) > smooth_win_size:
                smoothed_frequencies = {}
                for symbol in DNA_ALPHABET:
                    smoothed_frequencies[symbol] = self.apply_smoothing(
                        frequencies[symbol], smooth_win_size
                    )
                
                offset = (smooth_win_size - 1) // 2
                adjusted_positions = window_positions[offset:-offset]
                
                self.plot_frequencies(
                    smoothed_frequencies, 
                    adjusted_positions, 
                    first_record.id, 
                    is_smoothed, 
                    smooth_win_size
                )
            else:
                self.plot_frequencies(
                    frequencies, 
                    window_positions, 
                    first_record.id
                )

        except StopIteration:
            messagebox.showerror("Error", "The selected file is empty or not a valid FASTA file.")
        except Exception as e:
            messagebox.showerror("An Error Occurred", f"An unexpected error occurred: {e}")
            self.setup_initial_chart()

    def sliding_window_analysis(self, sequence):
        frequencies = {symbol: [] for symbol in DNA_ALPHABET}
        window_positions = []
        for i in range(len(sequence) - WINDOW_SIZE + 1):
            window = sequence[i : i + WINDOW_SIZE]
            for symbol in DNA_ALPHABET:
                count = window.count(symbol)
                relative_freq = count / WINDOW_SIZE
                frequencies[symbol].append(relative_freq)
            window_positions.append(i)
        return frequencies, window_positions

    def plot_frequencies(self, frequencies, positions, sequence_id, is_smoothed=False, smoothing_window=0):
        self.ax.clear()
        colors = {'A': 'blue', 'C': 'green', 'G': 'orange', 'T': 'red'}
        
        for symbol in DNA_ALPHABET:
            self.ax.plot(positions, frequencies[symbol], label=f"Freq of '{symbol}'", color=colors[symbol])

        title = f"Nucleotide Frequencies for: {sequence_id}"
        if is_smoothed:
            title += f" (Smoothed, Window: {smoothing_window})"
        
        self.ax.set_title(title)
        self.ax.set_xlabel("Window Start Position")
        self.ax.set_ylabel("Relative Frequency")
        self.ax.set_ylim(0, 1)
        self.ax.legend()
        self.ax.grid(True, linestyle='--', alpha=0.6)
        
        self.canvas.draw()

def main():
    root = tk.Tk()
    app = FastaAnalyzerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()