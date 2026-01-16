# A square matrix of an arbitrary size and the corresponding initial vector are given. 
# Implement a software application that makes a prediction on a total of 5 discrete steps, using this matrix and the corresponding vector.
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class MatrixPredictionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Matrix Prediction System - 5 Steps")
        self.root.geometry("1200x800")
        
        self.matrix_size = tk.IntVar(value=3)
        self.predictions = []
        
        self.create_widgets()
        self.create_matrix_inputs()
        
    def create_widgets(self):
        left_frame = ttk.Frame(self.root, padding="10")
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        right_frame = ttk.Frame(self.root, padding="10")
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=2)
        self.root.rowconfigure(0, weight=1)
        
        size_frame = ttk.LabelFrame(left_frame, text="Matrix Size", padding="10")
        size_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(size_frame, text="Size (n×n):").grid(row=0, column=0, sticky=tk.W)
        size_spinbox = ttk.Spinbox(size_frame, from_=2, to=10, textvariable=self.matrix_size, 
                                   width=10, command=self.create_matrix_inputs)
        size_spinbox.grid(row=0, column=1, padx=5)
        
        ttk.Button(size_frame, text="Update Size", 
                  command=self.create_matrix_inputs).grid(row=0, column=2, padx=5)
        
        self.matrix_frame = ttk.LabelFrame(left_frame, text="Transition Matrix", padding="10")
        self.matrix_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        self.vector_frame = ttk.LabelFrame(left_frame, text="Initial Vector", padding="10")
        self.vector_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)
        
        button_frame = ttk.Frame(left_frame)
        button_frame.grid(row=3, column=0, sticky=(tk.W, tk.E), pady=10)
        
        ttk.Button(button_frame, text="Run Prediction", 
                  command=self.run_prediction).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", 
                  command=self.clear_all).pack(side=tk.LEFT, padx=5)
        
        results_frame = ttk.LabelFrame(left_frame, text="Results", padding="10")
        results_frame.grid(row=4, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        self.results_text = scrolledtext.ScrolledText(results_frame, height=15, width=50)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        
        left_frame.rowconfigure(1, weight=1)
        left_frame.rowconfigure(4, weight=1)
        
        plot_frame = ttk.LabelFrame(right_frame, text="Prediction Plot", padding="10")
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.ax.set_xlabel('Step')
        self.ax.set_ylabel('Value')
        self.ax.set_title('State Evolution Over 5 Discrete Steps')
        self.ax.grid(True, alpha=0.3)
        
    def create_matrix_inputs(self):
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()
        for widget in self.vector_frame.winfo_children():
            widget.destroy()
        
        n = self.matrix_size.get()
        
        self.matrix_entries = []
        for i in range(n):
            row_entries = []
            for j in range(n):
                entry = ttk.Entry(self.matrix_frame, width=8)
                entry.grid(row=i, column=j, padx=2, pady=2)
                entry.insert(0, "0.0")
                row_entries.append(entry)
            self.matrix_entries.append(row_entries)
        
        self.vector_entries = []
        for i in range(n):
            ttk.Label(self.vector_frame, text=f"v[{i}]:").grid(row=i, column=0, sticky=tk.W, pady=2)
            entry = ttk.Entry(self.vector_frame, width=12)
            entry.grid(row=i, column=1, padx=5, pady=2)
            entry.insert(0, "1.0" if i == 0 else "0.0")
            self.vector_entries.append(entry)
    
    def get_matrix(self):
        n = self.matrix_size.get()
        matrix = np.zeros((n, n))
        try:
            for i in range(n):
                for j in range(n):
                    matrix[i, j] = float(self.matrix_entries[i][j].get())
            return matrix
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers in the matrix.")
            return None
    
    def get_vector(self):
        n = self.matrix_size.get()
        vector = np.zeros(n)
        try:
            for i in range(n):
                vector[i] = float(self.vector_entries[i].get())
            return vector
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers in the vector.")
            return None
    
    def predict_steps(self, matrix, initial_vector, num_steps=5):
        predictions = [initial_vector]
        current_vector = initial_vector.copy()
        
        for step in range(num_steps):
            current_vector = matrix @ current_vector
            predictions.append(current_vector)
        
        return predictions
    
    def run_prediction(self):
        matrix = self.get_matrix()
        initial_vector = self.get_vector()
        
        if matrix is None or initial_vector is None:
            return
        
        self.predictions = self.predict_steps(matrix, initial_vector, num_steps=5)
        
        self.display_results()
        
        self.update_plot()
    
    def display_results(self):
        self.results_text.delete(1.0, tk.END)
        
        self.results_text.insert(tk.END, "="*50 + "\n")
        self.results_text.insert(tk.END, "PREDICTION RESULTS\n")
        self.results_text.insert(tk.END, "="*50 + "\n\n")
        
        for step, vector in enumerate(self.predictions):
            self.results_text.insert(tk.END, f"Step {step}:\n")
            self.results_text.insert(tk.END, f"  Vector: {np.array2string(vector, precision=4, suppress_small=True)}\n")
            if step > 0:
                self.results_text.insert(tk.END, f"  Norm: {np.linalg.norm(vector):.4f}\n")
            self.results_text.insert(tk.END, "\n")
        
        self.results_text.insert(tk.END, "="*50 + "\n")
    
    def update_plot(self):
        self.ax.clear()
        
        predictions_array = np.array(self.predictions)
        num_steps = len(self.predictions) - 1
        num_states = predictions_array.shape[1]
        
        colors = plt.cm.tab10(np.linspace(0, 1, num_states))
        for i in range(num_states):
            self.ax.plot(range(num_steps + 1), predictions_array[:, i], 
                        marker='o', label=f'State {i+1}', color=colors[i], linewidth=2)
        
        self.ax.set_xlabel('Step', fontsize=12)
        self.ax.set_ylabel('Value', fontsize=12)
        self.ax.set_title('State Evolution Over 5 Discrete Steps', fontsize=14, fontweight='bold')
        self.ax.legend(loc='best')
        self.ax.grid(True, alpha=0.3)
        
        self.canvas.draw()
    
    def clear_all(self):
        self.create_matrix_inputs()
        self.results_text.delete(1.0, tk.END)
        self.ax.clear()
        self.ax.set_xlabel('Step')
        self.ax.set_ylabel('Value')
        self.ax.set_title('State Evolution Over 5 Discrete Steps')
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw()
    
    def load_example(self):
        self.matrix_size.set(3)
        self.create_matrix_inputs()
        
        example_matrix = [
            [0.7, 0.2, 0.1],
            [0.1, 0.6, 0.3],
            [0.2, 0.2, 0.6]
        ]
        
        example_vector = [1.0, 0.0, 0.0]
        
        for i in range(3):
            for j in range(3):
                self.matrix_entries[i][j].delete(0, tk.END)
                self.matrix_entries[i][j].insert(0, str(example_matrix[i][j]))
        
        for i in range(3):
            self.vector_entries[i].delete(0, tk.END)
            self.vector_entries[i].insert(0, str(example_vector[i]))
        
        messagebox.showinfo("Example Loaded", "Example 3x3 transition matrix and initial vector loaded!")

def main():
    root = tk.Tk()
    app = MatrixPredictionGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()