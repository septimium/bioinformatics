import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib. backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LinearSegmentedColormap


class NeedlemanWunschAligner:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("DNA Sequence Aligner - Needleman-Wunsch Algorithm")
        self.root.geometry("1100x700")
        
        # Default values
        self.seq1 = tk.StringVar(value="ACCGTGAAGCCAATAC")
        self.seq2 = tk.StringVar(value="AGCGTGCAGCCAATAC")
        self.gap_penalty = tk. IntVar(value=0)
        self.match_score = tk.IntVar(value=1)
        self.mismatch_score = tk.IntVar(value=-1)
        
        # Options
        self.plot_traceback = tk.BooleanVar(value=False)
        self.plot_grid = tk.BooleanVar(value=True)
        self.show_diagonal = tk. BooleanVar(value=False)
        
        self. score_matrix = None
        self.traceback_matrix = None
        
        self.setup_ui()
        
    def setup_ui(self):
        # Main frame
        main_frame = ttk. Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk. W, tk.E, tk.N, tk.S))
        
        # Left panel
        left_panel = ttk.Frame(main_frame)
        left_panel.grid(row=0, column=0, sticky=(tk. W, tk.N), padx=5)
        
        # Sequences frame
        seq_frame = ttk.LabelFrame(left_panel, text="Sequences", padding="5")
        seq_frame.grid(row=0, column=0, sticky=(tk.W, tk. E), pady=5)
        
        ttk. Label(seq_frame, text="Sq 1 =").grid(row=0, column=0, sticky=tk. W)
        ttk.Entry(seq_frame, textvariable=self.seq1, width=25).grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(seq_frame, text="Sq 2 =").grid(row=1, column=0, sticky=tk.W)
        ttk. Entry(seq_frame, textvariable=self.seq2, width=25).grid(row=1, column=1, padx=5, pady=2)
        
        # Parameters and Options frame
        params_options_frame = ttk.Frame(left_panel)
        params_options_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Parameters frame
        params_frame = ttk.LabelFrame(params_options_frame, text="Parameters", padding="5")
        params_frame.grid(row=0, column=0, sticky=(tk.W, tk.N), padx=(0, 10))
        
        ttk.Label(params_frame, text="Gap =").grid(row=0, column=0, sticky=tk. W)
        ttk.Entry(params_frame, textvariable=self. gap_penalty, width=5).grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="Match =").grid(row=1, column=0, sticky=tk.W)
        ttk.Entry(params_frame, textvariable=self.match_score, width=5).grid(row=1, column=1, padx=5, pady=2)
        
        ttk.Label(params_frame, text="MMatch =").grid(row=2, column=0, sticky=tk.W)
        ttk.Entry(params_frame, textvariable=self.mismatch_score, width=5).grid(row=2, column=1, padx=5, pady=2)
        
        # Options frame
        options_frame = ttk. LabelFrame(params_options_frame, text="Options", padding="5")
        options_frame.grid(row=0, column=1, sticky=(tk.W, tk.N))
        
        ttk.Checkbutton(options_frame, text="Plot TraceBack", variable=self.plot_traceback).grid(row=0, column=0, sticky=tk.W)
        ttk. Checkbutton(options_frame, text="Plot grid", variable=self. plot_grid).grid(row=1, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Show diagonal", variable=self.show_diagonal).grid(row=2, column=0, sticky=tk.W)
        
        # Align button
        ttk.Button(left_panel, text="Align", command=self.align).grid(row=2, column=0, pady=10)
        
        # Presets frame
        presets_frame = ttk.LabelFrame(left_panel, text="Presets", padding="5")
        presets_frame.grid(row=3, column=0, sticky=(tk.W, tk.E), pady=5)
        
        for i in range(8):
            ttk.Button(presets_frame, text=f"Setting {i+1}", width=10).grid(row=i, column=0, padx=2, pady=1)
            ttk.Button(presets_frame, text=f"Setting {i+9}", width=10).grid(row=i, column=1, padx=2, pady=1)
        
        # Right panel - Visualizations
        right_panel = ttk. Frame(main_frame)
        right_panel.grid(row=0, column=1, sticky=(tk.W, tk. E, tk.N, tk.S), padx=10)
        
        # Matrix visualization frame
        viz_frame = ttk.Frame(right_panel)
        viz_frame.grid(row=0, column=0, columnspan=2)
        
        # Alignment matrix visualization
        matrix_frame = ttk.LabelFrame(viz_frame, text="Alignment matrix", padding="5")
        matrix_frame.grid(row=0, column=0, padx=5)
        
        self.fig1, self.ax1 = plt.subplots(figsize=(2.8, 2.2))
        self.fig1.tight_layout()
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=matrix_frame)
        self.canvas1.get_tk_widget().pack()
        
        # Traceback visualization
        traceback_frame = ttk.LabelFrame(viz_frame, text="Traceback path", padding="5")
        traceback_frame.grid(row=0, column=1, padx=5)
        
        self.fig2, self.ax2 = plt.subplots(figsize=(2.8, 2.2))
        self.fig2.tight_layout()
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=traceback_frame)
        self.canvas2.get_tk_widget().pack()
        
        # Results text
        results_frame = ttk.LabelFrame(right_panel, text="Show Alignment:", padding="5")
        results_frame. grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        
        self.results_text = tk.Text(results_frame, height=12, width=55, font=("Courier", 10))
        self.results_text.pack(side=tk.LEFT, fill=tk. BOTH, expand=True)
        
        scrollbar = ttk. Scrollbar(results_frame, orient=tk.VERTICAL, command=self.results_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.results_text.configure(yscrollcommand=scrollbar. set)
        
    def needleman_wunsch(self, seq1, seq2, gap, match, mismatch):
        n = len(seq1) + 1
        m = len(seq2) + 1
        
        # Initialize score matrix
        score_matrix = np.zeros((n, m))
        traceback_matrix = np.zeros((n, m), dtype=int)
        
        # Initialize first row and column
        for i in range(n):
            score_matrix[i][0] = i * gap
            traceback_matrix[i][0] = 1  # Up
        for j in range(m):
            score_matrix[0][j] = j * gap
            traceback_matrix[0][j] = 2  # Left
        
        traceback_matrix[0][0] = 0
        
        # Fill the matrix
        for i in range(1, n):
            for j in range(1, m):
                if seq1[i-1] == seq2[j-1]: 
                    diag_score = score_matrix[i-1][j-1] + match
                else:
                    diag_score = score_matrix[i-1][j-1] + mismatch
                
                up_score = score_matrix[i-1][j] + gap
                left_score = score_matrix[i][j-1] + gap
                
                score_matrix[i][j] = max(diag_score, up_score, left_score)
                
                if score_matrix[i][j] == diag_score:
                    traceback_matrix[i][j] = 0  # Diagonal
                elif score_matrix[i][j] == up_score:
                    traceback_matrix[i][j] = 1  # Up
                else:
                    traceback_matrix[i][j] = 2  # Left
        
        return score_matrix, traceback_matrix
    
    def traceback(self, seq1, seq2, traceback_matrix):
        aligned1 = ""
        aligned2 = ""
        match_string = ""
        
        i, j = len(seq1), len(seq2)
        path = [(i, j)]
        
        while i > 0 or j > 0:
            if traceback_matrix[i][j] == 0:  # Diagonal
                aligned1 = seq1[i-1] + aligned1
                aligned2 = seq2[j-1] + aligned2
                if seq1[i-1] == seq2[j-1]: 
                    match_string = "|" + match_string
                else:
                    match_string = " " + match_string
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 1:  # Up
                aligned1 = seq1[i-1] + aligned1
                aligned2 = "-" + aligned2
                match_string = " " + match_string
                i -= 1
            else:   # Left
                aligned1 = "-" + aligned1
                aligned2 = seq2[j-1] + aligned2
                match_string = " " + match_string
                j -= 1
            path.append((i, j))
        
        return aligned1, aligned2, match_string, path
    
    def align(self):
        seq1 = self.seq1.get().upper()
        seq2 = self. seq2.get().upper()
        gap = self.gap_penalty.get()
        match = self.match_score.get()
        mismatch = self.mismatch_score.get()
        
        # Perform alignment
        self.score_matrix, self.traceback_matrix = self.needleman_wunsch(
            seq1, seq2, gap, match, mismatch
        )
        
        aligned1, aligned2, match_string, path = self.traceback(
            seq1, seq2, self.traceback_matrix
        )
        
        # Calculate statistics
        matches = match_string.count("|")
        length = len(aligned1)
        similarity = (matches / length) * 100 if length > 0 else 0
        
        # Update visualizations
        self. plot_score_matrix()
        self.plot_traceback_path(path, len(seq1) + 1, len(seq2) + 1)
        
        # Update results text
        self.results_text. delete(1.0, tk.END)
        self.results_text.insert(tk.END, f"{aligned1}\n")
        self.results_text.insert(tk.END, f"{match_string}\n")
        self.results_text.insert(tk.END, f"{aligned2}\n\n")
        self.results_text. insert(tk.END, f"Matches = {matches}\n")
        self.results_text.insert(tk.END, f"Length = {length}\n\n")
        self.results_text. insert(tk.END, f"Similarity = {similarity:.0f} %\n\n")
        self.results_text. insert(tk.END, f"Tracing back:  M[{len(seq1)},{len(seq2)}]\n\n")
        
        # Display score matrix
        self. results_text.insert(tk.END, "Score Matrix:\n")
        header = "      " + "  ".join([f"{c:>3}" for c in " " + seq2]) + "\n"
        self.results_text.insert(tk.END, header)
        for i, row in enumerate(self.score_matrix):
            if i == 0:
                row_label = " "
            else: 
                row_label = seq1[i-1]
            row_str = f"{row_label: >3} " + "  ".join([f"{int(v):>3}" for v in row]) + "\n"
            self.results_text.insert(tk.END, row_str)
    
    def plot_score_matrix(self):
        self.ax1.clear()
        
        # Create custom colormap (dark blue to red)
        colors = ['#000033', '#330066', '#660066', '#990033', '#CC0000', '#FF3333']
        cmap = LinearSegmentedColormap.from_list('custom', colors)
        
        self.ax1.imshow(self.score_matrix, cmap=cmap, aspect='auto')
        
        if self. plot_grid.get():
            self.ax1.set_xticks(np.arange(-0.5, self. score_matrix.shape[1], 1), minor=True)
            self.ax1.set_yticks(np.arange(-0.5, self.score_matrix.shape[0], 1), minor=True)
            self.ax1.grid(which='minor', color='white', linestyle='-', linewidth=0.5, alpha=0.3)
        
        self.ax1.set_xticks([])
        self.ax1.set_yticks([])
        self.fig1.tight_layout()
        self.canvas1.draw()
    
    def plot_traceback_path(self, path, n, m):
        self.ax2.clear()
        
        # Create path matrix
        path_matrix = np. zeros((n, m))
        for (i, j) in path:
            path_matrix[i][j] = 1
        
        # Create custom colormap
        cmap = plt.cm.colors. ListedColormap(['#FFFFCC', '#CC0000'])
        
        self.ax2.imshow(path_matrix, cmap=cmap, aspect='auto')
        
        # Add grid
        self.ax2.set_xticks(np.arange(-0.5, m, 1), minor=True)
        self.ax2.set_yticks(np.arange(-0.5, n, 1), minor=True)
        self.ax2.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
        
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])
        self.fig2.tight_layout()
        self.canvas2.draw()
    
    def run(self):
        self.root.mainloop()


if __name__ == "__main__": 
    app = NeedlemanWunschAligner()
    app.run()