import os
import sys
import ctypes
import tkinter as tk
import ttkbootstrap as tb
from ttkbootstrap.constants import *
from tkinter import ttk, messagebox, filedialog
import re
import threading
import time
import numpy as np
from pathlib import Path
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure



class CCPApp(tb.Window):  # Bootstrap themed window
    def __init__(self):
        super().__init__(themename="cyborg")  

        self.title("Plasma Reactor Simulation")
        self.update_idletasks()
        width = self.winfo_screenwidth()
        height = self.winfo_screenheight()
        self.geometry(f"{min(1000, width - 100)}x{min(700, height - 100)}")
        self.resizable(True, True)
        self.enable_Optimisation = tk.BooleanVar(value=True)

        # Use modern font
        self.option_add("*Font", ("Segoe UI", 10))

        self.param_vars = {
            'power': tk.DoubleVar(value=500.0),
            'pressure': tk.DoubleVar(value=10.0),
            'gap_length': tk.DoubleVar(value=0.025),
            'temperature': tk.DoubleVar(value=350.0),
            'turns': tk.DoubleVar(value=5),  # default to 5 turn coil
            'target_efficiency' : tk.DoubleVar(value=70),
        }
        self.param_min_vars = {}
        self.param_max_vars = {}


        # Load the shared library
        self.project_dir = Path(__file__).parent.resolve()
        self.ccp = self.load_ccp_library()

        self.output_folder = ""  

        #simulation locking for threading
        self.simulation_lock = threading.Lock()
        self.simulation_running = False  

        self.plot_index = 0

        # GUI
        self.create_widgets()

    def select_output_folder(self):
        folder = filedialog.askdirectory(title="Select Output Folder")
        if folder:
            self.output_folder = folder
            self.output_text.insert("end", f">> Selected working folder: {folder}\n")


    def get_latest_results_dir(self):
        if not self.output_folder:
            messagebox.showerror("Data Error", "Output folder not set.")
            return self.project_dir

        output_path = Path(self.output_folder)
        pattern = re.compile(r"run_(\d+)")
        runs = []

        for d in output_path.iterdir():
            if d.is_dir():
                match = pattern.match(d.name)
                if match:
                    runs.append((int(match.group(1)), d))

        if runs:
            latest_run = max(runs, key=lambda x: x[0])[1]
            return latest_run
        else:
            messagebox.showerror("Data Error", f"No run folders found in: {output_path}")
            return output_path



    def load_ccp_library(self):
        base_dir = Path(__file__).parent.resolve()
        libname = "simVI.dll" if os.name == "nt" else "libccpsim.so"
        libpath = base_dir / libname

        if not libpath.exists():
            messagebox.showerror("Error", f"Simulation library not found:\n{libpath}")
            sys.exit(1)

        try:
            ccp = ctypes.CDLL(str(libpath))
            ccp.run_simulation.argtypes = [
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char,
                ctypes.c_int,
                ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                ctypes.c_double, ctypes.c_double
            ]

            ccp.run_simulation.restype = ctypes.c_int
            return ccp
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load library:\n{e}")
            sys.exit(1)


    def create_widgets(self):
        main_frame = tb.Frame(self, padding=10)
        main_frame.pack(fill="both", expand=True)

        # Notebook for tabs
        self.notebook = tb.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True)
        self.notebook.enable_traversal()

        control_tab = tb.Frame(self.notebook)
        self.notebook.add(control_tab, text="Controls")
        control_tab.columnconfigure(0, weight=1)
        control_tab.columnconfigure(1, weight=1)
        control_tab.rowconfigure(3, weight=1)

        # Reactor Parameters Frame
        param_frame = tb.Labelframe(control_tab, text="Reactor Parameters", padding=10)
        param_frame.grid(row=0, column=0, columnspan=2, sticky="ew", pady=5)
        param_frame.columnconfigure(1, weight=1)

        # Add column labels above value/min/max
        tb.Label(param_frame, text="Value").grid(row=0, column=1, sticky='s')
        tb.Label(param_frame, text="Min").grid(row=0, column=2, sticky='s')
        tb.Label(param_frame, text="Max").grid(row=0, column=3, sticky='s')

        params = [
            ('Power (W):', 'power'),
            ('Pressure (Pa):', 'pressure'),
            ('Gap Length (m):', 'gap_length'),
        ]
        default_bounds = {
            'power':       (5.0, 5000.0),
            'pressure':    (3.0, 15.0),
            'gap_length':  (0.025, 0.1),
        }
        # Add Number of Coil Turns and Gas Temp(no min/max)
        row_offset = len(params) + 1
        # Temperature input
        tb.Label(param_frame, text="Background Gas Temperature (K):").grid(row=row_offset, column=0, sticky='w', pady=6, padx=5)
        tb.Entry(param_frame, textvariable=self.param_vars['temperature'], width=10).grid(row=row_offset, column=1, padx=5)

        # Number of coil turns
        tb.Label(param_frame, text="Number of Coils:").grid(row=row_offset+1, column=0, sticky='w', pady=6, padx=5)
        tb.Entry(param_frame, textvariable=self.param_vars['turns'], width=10).grid(row=row_offset+1, column=1, padx=5)

        for i, (label, var_name) in enumerate(params, start=1):
            tb.Label(param_frame, text=label).grid(row=i, column=0, sticky='w', pady=6, padx=5)

            # Main value
            tb.Entry(param_frame, textvariable=self.param_vars[var_name], width=10).grid(row=i, column=1, padx=5)

            # Min value
            min_default, max_default = default_bounds[var_name]
            min_var = tk.DoubleVar(value=min_default)
            self.param_min_vars[var_name] = min_var
            tb.Entry(param_frame, textvariable=min_var, width=10).grid(row=i, column=2, padx=5)

            # Max value
            max_var = tk.DoubleVar(value=max_default)
            self.param_max_vars[var_name] = max_var
            tb.Entry(param_frame, textvariable=max_var, width=10).grid(row=i, column=3, padx=5)

        # Plasma Type Frame
        plasma_frame = tb.Labelframe(control_tab, text="Plasma Type", padding=10)
        plasma_frame.grid(row=1, column=0, sticky="w", pady=5, padx=5)

        self.sim_type_var = tk.StringVar(value='C')
        tb.Radiobutton(plasma_frame, text="CCP", variable=self.sim_type_var, value='C').pack(side="left", padx=5)
        tb.Radiobutton(plasma_frame, text="ICP", variable=self.sim_type_var, value='I').pack(side="left", padx=5)

        def toggle_opt_dropdown_state(*args):
            state = 'readonly' if self.enable_Optimisation.get() else 'disabled'
            self.opt_type_dropdown.configure(state=state)

        self.enable_Optimisation.trace_add('write', lambda *args: toggle_opt_dropdown_state())

        # Optimisation Parameter Frame
        opt_frame = tb.Labelframe(control_tab, text="Optimisation Parameter", padding=10)
        opt_frame.grid(row=1, column=1, sticky="ew", pady=5, padx=5)
        opt_frame.columnconfigure(0, weight=1)

        self.opt_type_var = tk.StringVar()
        self.opt_type_dropdown = tb.Combobox(
            opt_frame,
            textvariable=self.opt_type_var,
            values=["Power", "Pressure", "Gap Length"],
            state="readonly",
        )
        self.opt_type_dropdown.current(0)
        self.opt_type_dropdown.pack(fill="x")

        # optimisation check box
        self.opt_checkbox = tb.Checkbutton(
            opt_frame,
            text="Enable Optimisation",
            variable=self.enable_Optimisation,
            onvalue=True,
            offvalue=False
        )
        self.opt_checkbox.pack(anchor="w", pady=5)


        # Row for Cycles and Target Efficiency side-by-side
        cycles_eff_frame = tb.Frame(control_tab)
        cycles_eff_frame.grid(row=2, column=0, columnspan=2, sticky="w", padx=5, pady=10)

        # Cycles input
        tb.Label(cycles_eff_frame, text="Cycles:").grid(row=0, column=0, sticky="w")
        self.cycles_var = tk.StringVar(value="1000")
        tb.Entry(cycles_eff_frame, textvariable=self.cycles_var, width=10).grid(row=0, column=1, padx=(5, 15))

        # Target Efficiency input
        tb.Label(cycles_eff_frame, text="Target Efficiency (%):").grid(row=0, column=2, sticky="w")
        tb.Entry(cycles_eff_frame, textvariable=self.param_vars['target_efficiency'], width=10).grid(row=0, column=3, padx=(5, 0))

        # Buttons Frame
        btn_frame = tb.Frame(control_tab)
        btn_frame.grid(row=4, column=0, columnspan=2, pady=10, sticky="ew")
        btn_frame.columnconfigure(5, weight=1, uniform="a")

        tb.Button(btn_frame, text="Initialise", command=self.run_initialisation).grid(row=0, column=0, padx=6, sticky="ew")
        tb.Button(btn_frame, text="Run Simulation", command=self.run_simulation, style="Accent.TButton").grid(row=0, column=1, padx=6, sticky="ew")
        tb.Button(btn_frame, text="Visualise Results", command=self.show_visualisations).grid(row=0, column=2, padx=6, sticky="ew")
        tb.Button(btn_frame, text="Clear Log", command=self.clear_log).grid(row=0, column=3, padx=6, sticky="ew")
        btn_frame.columnconfigure(4, weight=1, uniform="a")
        tb.Button(btn_frame, text="Reset Defaults", command=self.reset_to_defaults).grid(
            row=0, column=4, padx=6, sticky="ew"
        )
        tb.Button(btn_frame, text="Select Folder", command=self.select_output_folder).grid(row=0, column=5, padx=6, sticky="ew")



        # Progress bar
        self.progress = tb.Progressbar(control_tab, mode="indeterminate", bootstyle="info-striped", length=200)
        self.progress.grid(row=5, column=0, columnspan=2, pady=(0, 10))

        # Output Frame with Scrollbar
        output_frame = tb.Labelframe(control_tab, text="Simulation Output", padding=10)
        output_frame.grid(row=3, column=0, columnspan=2, sticky="nsew", pady=5)
        control_tab.rowconfigure(3, weight=1)

        self.output_text = tb.ScrolledText(output_frame, height=12, wrap="word", font=("Segoe UI", 10))
        self.output_text.pack(side="left", fill="both", expand=True)

        scrollbar = tb.Scrollbar(output_frame, command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.output_text.configure(yscrollcommand=scrollbar.set)   
    
    def start_progress(self):
        self.progress.start(10)  # Speed (ms per step)

    def stop_progress(self):
        self.progress.stop()

        
    def reset_to_defaults(self):
        default_bounds = {
            'power':       (5.0, 5000.0),
            'pressure':    (3.0, 15.0),
            'gap_length':  (0.025, 0.1)
        }
        
        # Default values for parameters
        default_values = {
            'power': 500.0,
            'pressure': 10.0,
            'gap_length': 0.025,
            'temperature': 350.0,
            'turns': 5.0,
            'target_efficiency': 70.0
        }

        # Reset all parameter values
        for key in self.param_vars:
            if key in default_values:
                self.param_vars[key].set(default_values[key])
            
            # Reset min/max bounds (only for parameters that have bounds)
            if key in default_bounds:
                default_min, default_max = default_bounds[key]
                self.param_min_vars[key].set(default_min)
                self.param_max_vars[key].set(default_max)
        
        # Reset other GUI elements to defaults
        self.sim_type_var.set('C')  # Default to CCP
        self.cycles_var.set("1000")  # Default cycles
        self.param_vars['target_efficiency'].set(70.0)
        self.enable_Optimisation.set(True)  # Default Optimisation enabled
        self.opt_type_var.set("Power")  # Default Optimisation parameter
        
        # Clear output log
        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, ">> Parameters reset to defaults\n")

    def create_parameter_entry(self, parent, var, row, column):
        entry = tb.Entry(parent, textvariable=var, width=15)
        entry.grid(row=row, column=column, sticky="ew", padx=5, pady=3)
        parent.columnconfigure(column, weight=1)
        return entry

    def clear_log(self):
        self.output_text.delete("1.0", tk.END) 
        
    def embed_responsive_plot(self, fig, parent, filename_base, index):
        fig.set_size_inches(6, 3.5)
        plot_frame = tb.Frame(parent)
        plot_frame.pack(fill="both", expand=True, padx=10, pady=10)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        widget = canvas.get_tk_widget()
        widget.grid(row=0, column=0, sticky="nsew")

        canvas.draw()

        # Export buttons
        btn_frame = tb.Frame(plot_frame)
        btn_frame.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        btn_frame.columnconfigure(0, weight=1)

        def export(fmt):
            filetypes = [(f"{fmt.upper()} files", f"*.{fmt}"), ("All files", "*.*")]
            default_name = f"{filename_base}.{fmt}"
            filepath = filedialog.asksaveasfilename(
                defaultextension=f".{fmt}",
                filetypes=filetypes,
                initialfile=default_name,
                title=f"Save as {fmt.upper()}"
            )
            if filepath:
                try:
                    fig.savefig(filepath, bbox_inches='tight')
                    messagebox.showinfo("Export", f"Saved: {os.path.basename(filepath)}")
                except Exception as e:
                    messagebox.showerror("Export Error", f"Failed to save file:\n{str(e)}")

        tb.Button(btn_frame, text="Export PNG", command=lambda: export("png")).pack(side='left', padx=5)
        tb.Button(btn_frame, text="Export PDF", command=lambda: export("pdf")).pack(side='left')

    def show_visualisations(self):
        # Find existing Results tab
        results_tab_id = None
        for tab_id in range(self.notebook.index("end")):
            if self.notebook.tab(tab_id, "text") == "Results":
                results_tab_id = tab_id
                break
        
        if results_tab_id is not None:
            # Remove existing Results tab
            self.notebook.forget(results_tab_id)
        
        # Create new Results tab
        results_tab = tb.Frame(self.notebook)
        results_tab.rowconfigure(0, weight=1)   
        results_tab.columnconfigure(0, weight=1)  
        self.notebook.add(results_tab, text="Results")
        self.create_visualisations(results_tab)
        
        # Select the new Results tab
        for tab_id in range(self.notebook.index("end")):
            if self.notebook.tab(tab_id, "text") == "Results":
                self.notebook.select(tab_id)
                break


    def tab_exists(self, tab_text):
        for tab_id in range(self.notebook.index("end")):
            if self.notebook.tab(tab_id, "text") == tab_text:
                return True
        return False

    def plot_heating_efficiency(self, parent):
        """Plot heating efficiency across cycles for current run and optionally compare multiple runs"""
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            # Get current run data
            current_entries = self.read_heating_efficiency_log(self.data_dir)
            
            if current_entries:
                cycles = [entry['cycle'] for entry in current_entries]
                efficiencies = [entry['efficiency'] for entry in current_entries]
                
                # Check if this is Optimisation mode (has parameter data)
                has_params = 'param_name' in current_entries[0] and 'parameter' in current_entries[0]
                
                if has_params:
                    # Create dual y-axis plot for Optimisation runs
                    ax2 = ax.twinx()
                    
                    # Plot efficiency on left axis
                    line1 = ax.plot(cycles, efficiencies, 'b-o', label='Efficiency', markersize=4)
                    ax.set_ylabel('Efficiency (%)', color='b')
                    ax.tick_params(axis='y', labelcolor='b')
                    
                    # Plot parameter on right axis
                    param_name = current_entries[0]['param_name']
                    param_values = [entry['parameter'] for entry in current_entries]
                    line2 = ax2.plot(cycles, param_values, 'r-s', label=param_name, markersize=4)
                    ax2.set_ylabel(f'{param_name}', color='r')
                    ax2.tick_params(axis='y', labelcolor='r')
                    
                    # Combined legend
                    lines = line1 + line2
                    labels = [l.get_label() for l in lines]
                    ax.legend(lines, labels, loc='upper left')
                    
                    ax.set_title(f'Heating Efficiency Optimisation - {self.data_dir.name}')
                else:
                    # Simple efficiency plot
                    ax.plot(cycles, efficiencies, 'b-o', label='Efficiency', markersize=4)
                    ax.set_ylabel('Efficiency (%)')
                    ax.legend()
                    ax.set_title(f'Heating Efficiency - {self.data_dir.name}')
                
                ax.set_xlabel('Cycle')
                ax.grid(True, alpha=0.3)
                
                # Add final efficiency as text annotation
                final_eff = efficiencies[-1]
                ax.annotate(f'Final: {final_eff:.2f}%', 
                        xy=(cycles[-1], final_eff),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

            self.embed_responsive_plot(fig, frame, "heating_efficiency", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Heating efficiency plot failed: {str(e)}")

    def plot_efficiency_comparison(self, parent):
        """Compare heating efficiency across multiple runs"""
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            if not self.output_folder:
                ax.text(0.5, 0.5, 'No output folder selected', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_comparison", self.plot_index)
                return

            # Find all run folders
            output_path = Path(self.output_folder)
            pattern = re.compile(r"run_(\d+)")
            runs = []

            for d in output_path.iterdir():
                if d.is_dir():
                    match = pattern.match(d.name)
                    if match:
                        runs.append((int(match.group(1)), d))

            runs.sort(key=lambda x: x[0])  # Sort by run number

            if not runs:
                ax.text(0.5, 0.5, 'No run folders found', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_comparison", self.plot_index)
                return

            # Collect data from all runs
            run_numbers = []
            final_efficiencies = []
            Optimisation_params = []
            param_names = []

            for run_num, run_path in runs:
                entries = self.read_heating_efficiency_log(run_path)
                if entries:
                    run_numbers.append(run_num)
                    final_efficiencies.append(entries[-1]['efficiency'])
                    
                    # Check if this run has Optimisation parameters
                    if 'param_name' in entries[-1] and 'parameter' in entries[-1]:
                        Optimisation_params.append(entries[-1]['parameter'])
                        param_names.append(entries[-1]['param_name'])
                    else:
                        Optimisation_params.append(None)
                        param_names.append(None)

            if not run_numbers:
                ax.text(0.5, 0.5, 'No valid efficiency data found', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_comparison", self.plot_index)
                return

            # Check if we have Optimisation data
            has_opt_data = any(param is not None for param in Optimisation_params)
            
            if has_opt_data:
                # Create scatter plot with parameter values as colors
                unique_param_name = next(name for name in param_names if name is not None)
                valid_indices = [i for i, param in enumerate(Optimisation_params) if param is not None]
                
                if valid_indices:
                    valid_runs = [run_numbers[i] for i in valid_indices]
                    valid_effs = [final_efficiencies[i] for i in valid_indices]
                    valid_params = [Optimisation_params[i] for i in valid_indices]
                    
                    scatter = ax.scatter(valid_runs, valid_effs, c=valid_params, 
                                    cmap='viridis', s=60, alpha=0.7)
                    
                    # Add colorbar
                    cbar = fig.colorbar(scatter, ax=ax)
                    cbar.set_label(unique_param_name)
                    
                    ax.set_title(f'Heating Efficiency vs {unique_param_name} Across Runs')
                else:
                    ax.plot(run_numbers, final_efficiencies, 'o-', markersize=6)
                    ax.set_title('Heating Efficiency Across Runs')
            else:
                # Simple line plot
                ax.plot(run_numbers, final_efficiencies, 'o-', markersize=6)
                ax.set_title('Heating Efficiency Across Runs')

            ax.set_xlabel('Run Number')
            ax.set_ylabel('Final Efficiency (%)')
            ax.grid(True, alpha=0.3)

            # Add best efficiency annotation
            if final_efficiencies:
                best_idx = final_efficiencies.index(max(final_efficiencies))
                best_run = run_numbers[best_idx]
                best_eff = final_efficiencies[best_idx]
                
                ax.annotate(f'Best: Run {best_run}\n{best_eff:.2f}%', 
                        xy=(best_run, best_eff),
                        xytext=(10, 10), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7),
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

            self.embed_responsive_plot(fig, frame, "efficiency_comparison", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Efficiency comparison plot failed: {str(e)}")

    def plot_efficiency_statistics(self, parent):
        """Show statistical summary of heating efficiency across runs"""
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            if not self.output_folder:
                ax.text(0.5, 0.5, 'No output folder selected', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_statistics", self.plot_index)
                return

            # Collect all efficiency data
            output_path = Path(self.output_folder)
            pattern = re.compile(r"run_(\d+)")
            all_efficiencies = []

            for d in output_path.iterdir():
                if d.is_dir() and pattern.match(d.name):
                    entries = self.read_heating_efficiency_log(d)
                    if entries:
                        all_efficiencies.append(entries[-1]['efficiency'])

            if not all_efficiencies:
                ax.text(0.5, 0.5, 'No efficiency data found', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_statistics", self.plot_index)
                return

            # Create histogram
            ax.hist(all_efficiencies, bins=min(20, len(all_efficiencies)), 
                    alpha=0.7, color='skyblue', edgecolor='black')
            
            # Add statistics
            mean_eff = np.mean(all_efficiencies)
            std_eff = np.std(all_efficiencies)
            median_eff = np.median(all_efficiencies)
            
            ax.axvline(mean_eff, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_eff:.2f}%')
            ax.axvline(median_eff, color='green', linestyle='--', linewidth=2, label=f'Median: {median_eff:.2f}%')
            
            ax.set_xlabel('Heating Efficiency (%)')
            ax.set_ylabel('Number of Runs')
            ax.set_title(f'Distribution of Heating Efficiency\n(σ = {std_eff:.2f}%, n = {len(all_efficiencies)})')
            ax.legend()
            ax.grid(True, alpha=0.3)

            self.embed_responsive_plot(fig, frame, "efficiency_statistics", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Efficiency statistics plot failed: {str(e)}")

    def plot_efficiency_vs_parameter(self, parent):
        """Plot efficiency vs Optimisation parameter to show the relationship"""
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            if not self.output_folder:
                ax.text(0.5, 0.5, 'No output folder selected', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_vs_parameter", self.plot_index)
                return

            # Collect data from all runs that have Optimisation parameters
            output_path = Path(self.output_folder)
            pattern = re.compile(r"run_(\d+)")
            
            param_values = []
            efficiencies = []
            run_numbers = []
            param_name = None

            for d in output_path.iterdir():
                if d.is_dir():
                    match = pattern.match(d.name)
                    if match:
                        run_num = int(match.group(1))
                        entries = self.read_heating_efficiency_log(d)
                        if entries and 'param_name' in entries[-1] and 'parameter' in entries[-1]:
                            param_values.append(entries[-1]['parameter'])
                            efficiencies.append(entries[-1]['efficiency'])
                            run_numbers.append(run_num)
                            if param_name is None:
                                param_name = entries[-1]['param_name']

            if not param_values:
                ax.text(0.5, 0.5, 'No Optimisation parameter data found\n(Only available for Optimisation runs)', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "efficiency_vs_parameter", self.plot_index)
                return

            # Sort by parameter value for better visualization
            sorted_data = sorted(zip(param_values, efficiencies, run_numbers))
            param_values, efficiencies, run_numbers = zip(*sorted_data)

            # Create scatter plot with run numbers as colors
            scatter = ax.scatter(param_values, efficiencies, c=run_numbers, 
                            cmap='viridis', s=60, alpha=0.7, edgecolors='black')
            
            # Add colorbar for run numbers
            cbar = fig.colorbar(scatter, ax=ax)
            cbar.set_label('Run Number')

            # Add trend line if we have enough points
            if len(param_values) > 3:
                # Fit polynomial trend line
                z = np.polyfit(param_values, efficiencies, 2)  # 2nd degree polynomial
                p = np.poly1d(z)
                
                # Generate smooth curve
                x_smooth = np.linspace(min(param_values), max(param_values), 100)
                y_smooth = p(x_smooth)
                
                ax.plot(x_smooth, y_smooth, 'r--', alpha=0.8, linewidth=2, label='Trend')
                ax.legend()

            # Find and highlight optimal point
            best_idx = efficiencies.index(max(efficiencies))
            best_param = param_values[best_idx]
            best_eff = efficiencies[best_idx]
            best_run = run_numbers[best_idx]
            
            ax.annotate(f'Optimal: Run {best_run}\n{param_name} = {best_param:.3g}\nEff = {best_eff:.2f}%', 
                    xy=(best_param, best_eff),
                    xytext=(20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))

            ax.set_xlabel(f'{param_name}')
            ax.set_ylabel('Heating Efficiency (%)')
            ax.set_title(f'Heating Efficiency vs {param_name}')
            ax.grid(True, alpha=0.3)

            # Add some statistics to the plot
            correlation = np.corrcoef(param_values, efficiencies)[0, 1]
            ax.text(0.02, 0.98, f'Correlation: {correlation:.3f}', 
                    transform=ax.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

            self.embed_responsive_plot(fig, frame, "efficiency_vs_parameter", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Efficiency vs parameter plot failed: {str(e)}")

    def plot_parameter_evolution(self, parent):
        """Show how parameter values evolve during Optimisation within a single run"""
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            # Get current run data
            current_entries = self.read_heating_efficiency_log(self.data_dir)
            
            if not current_entries or 'param_name' not in current_entries[0]:
                ax.text(0.5, 0.5, 'No Optimisation parameter data found\n(Only available for Optimisation runs)', 
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
                self.embed_responsive_plot(fig, frame, "parameter_evolution", self.plot_index)
                return

            cycles = [entry['cycle'] for entry in current_entries]
            efficiencies = [entry['efficiency'] for entry in current_entries]
            param_values = [entry['parameter'] for entry in current_entries]
            param_name = current_entries[0]['param_name']

            # Create scatter plot with efficiency as color
            scatter = ax.scatter(cycles, param_values, c=efficiencies, 
                            cmap='RdYlGn', s=50, alpha=0.7, edgecolors='black')
            
            # Add colorbar for efficiency
            cbar = fig.colorbar(scatter, ax=ax)
            cbar.set_label('Efficiency (%)')

            # Add connecting line to show evolution
            ax.plot(cycles, param_values, 'b-', alpha=0.3, linewidth=1)

            # Highlight start and end points
            ax.scatter(cycles[0], param_values[0], color='red', s=100, marker='s', 
                    label=f'Start: {param_values[0]:.3g}', zorder=5)
            ax.scatter(cycles[-1], param_values[-1], color='green', s=100, marker='*', 
                    label=f'End: {param_values[-1]:.3g}', zorder=5)

            ax.set_xlabel('Cycle')
            ax.set_ylabel(f'{param_name}')
            ax.set_title(f'{param_name} Evolution During Optimisation - {self.data_dir.name}')
            ax.legend()
            ax.grid(True, alpha=0.3)

            # Add improvement annotation
            efficiency_improvement = efficiencies[-1] - efficiencies[0]
            param_change = param_values[-1] - param_values[0]
            
            ax.text(0.02, 0.98, 
                    f'Efficiency: {efficiencies[0]:.2f}% → {efficiencies[-1]:.2f}% ({efficiency_improvement:+.2f}%)\n'
                    f'{param_name}: {param_values[0]:.3g} → {param_values[-1]:.3g} ({param_change:+.3g})', 
                    transform=ax.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

            self.embed_responsive_plot(fig, frame, "parameter_evolution", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Parameter evolution plot failed: {str(e)}")

    # Modified create_visualisations method to include efficiency plots
    def create_visualisations(self, parent):
        self.data_dir = self.get_latest_results_dir()
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)

        container = tb.Frame(parent)
        container.grid(row=0, column=0, sticky="nsew")
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)

        canvas = tk.Canvas(container)
        scrollbar = tb.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")

        container.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        scrollable_frame = tb.Frame(canvas)
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        scrollable_frame.grid_columnconfigure(0, weight=1)
        scrollable_frame.grid_columnconfigure(1, weight=1)

        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

        def on_canvas_resize(event):
            canvas.itemconfig(canvas_window, width=event.width)

        canvas.bind("<Configure>", on_canvas_resize)
        scrollable_frame.columnconfigure(0, weight=1)

        # Add all plots including efficiency plots
        self.plot_index = 0
        self.plot_convergence(scrollable_frame)
        self.plot_index += 1
        self.plot_distributions(scrollable_frame)
        self.plot_index += 1
        self.plot_temp_avg_density(scrollable_frame)
        self.plot_index += 1
        self.plot_heating_efficiency(scrollable_frame)  # NEW: Single run efficiency
        self.plot_index += 1
        self.plot_efficiency_comparison(scrollable_frame)  # NEW: Multi-run comparison
        self.plot_index += 1
        self.plot_efficiency_vs_parameter(scrollable_frame)  # NEW: Efficiency vs parameter
        self.plot_index += 1
        self.plot_parameter_evolution(scrollable_frame)  # NEW: Parameter evolution in single run
        self.plot_index += 1
        self.plot_efficiency_statistics(scrollable_frame)  # NEW: Statistical summary


    def plot_convergence(self, parent):
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=150)
            ax = fig.add_subplot(111)

            data = np.loadtxt(self.data_dir / "conv.dat")
            ax.plot(data[:, 0], data[:, 1], label="Electrons", marker='o')
            ax.plot(data[:, 0], data[:, 2], label="Ions", marker='s')
            ax.set_xlabel("RF Cycles")
            ax.set_ylabel("Superparticles")
            ax.legend()
            ax.grid(True)

            self.embed_responsive_plot(fig, frame, "convergence", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Convergence plot failed: {str(e)}")


    def plot_distributions(self, parent):
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            data = np.loadtxt(self.data_dir / "ne_xt.dat")
            position = data[:, 0]
            time_steps = np.arange(1, data.shape[1])

            ax.imshow(data[:, 1:], aspect='auto',
                    extent=[time_steps[0], time_steps[-1], position[0], position[-1]],
                    origin='lower', cmap='plasma')
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Position (mm)")

            self.embed_responsive_plot(fig, frame, "distribution", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Distribution plot failed: {str(e)}")

    def plot_temp_avg_density(self, parent):
        try:
            frame = tb.Frame(parent)
            frame.pack(fill='both', expand=True, padx=5, pady=5)

            fig = Figure(dpi=100)
            ax = fig.add_subplot(111)

            data = np.loadtxt(self.data_dir / "density.dat")
            position = data[:, 0]
            ni_values = data[:, 1::2]
            ne_values = data[:, 2::2]

            ax.plot(position, np.mean(ne_values, axis=1), label="$\\langle n_e \\rangle$")
            ax.plot(position, np.mean(ni_values, axis=1), label="$\\langle n_i \\rangle$")
            ax.set_xlabel("Position between electrodes")
            ax.set_ylabel("Densities $n_e$, $n_i$")
            #ax.set_ylim(0,1e17)
            ax.set_title("Temporally Averaged Particle Density")
            ax.legend()
            ax.grid(True)

            self.embed_responsive_plot(fig, frame, "densities", self.plot_index)

        except Exception as e:
            messagebox.showerror("Plot Error", f"Temp. Avg. Density plot failed: {str(e)}")

    def get_simulation_params(self):
        try:
            sim_type = self.sim_type_var.get()  # 'C' = CCP, 'I' = ICP
            
            return {
                'power' : self.param_vars['power'].get(),
                'pressure': self.param_vars['pressure'].get(),
                'gap_length': self.param_vars['gap_length'].get(),
                'temperature': self.param_vars['temperature'].get(),
                'turns': self.param_vars['turns'].get(),
                'target_efficiency': self.param_vars['target_efficiency'].get()
            }
        except tk.TclError:
            messagebox.showerror("Input Error", "Invalid numeric values")
            return 
        
    def get_all_parameter_bounds(self):
        """Return dict of {parameter: {'value': ..., 'min': ..., 'max': ...}}"""
        bounds = {}
        try:
            for key in ['power', 'pressure', 'gap_length']:
                bounds[key] = {
                    "value": self.param_vars[key].get(),
                    "min": self.param_min_vars[key].get(),
                    "max": self.param_max_vars[key].get()
                }
        except tk.TclError as e:
            messagebox.showerror("Input Error", f"Invalid input in parameter fields:\n{e}")
            return None
        return bounds

    def run_simulation(self):
        if self.simulation_running:
            messagebox.showinfo("Simulation Busy", "A simulation is already running.")
            return

        # Start thread and pass opt_param
        threading.Thread(
            target=self._run_simulation_phase,
            args=(False, self.opt_type_var.get() if self.enable_Optimisation.get() else "None"),
            daemon=True
        ).start()

    def run_initialisation(self):
        if self.simulation_running:
            messagebox.showinfo("Simulation Busy", "A simulation is already running.")
            return
        threading.Thread(target=self._run_simulation_phase, args=(True,), daemon=True).start()


    def _run_simulation_phase(self, init=False, opt_param="None"):
        if not self.simulation_lock.acquire(blocking=False):
            print("Simulation already in progress.")
            return

        self.simulation_running = True
        self.start_progress()

        try:
            self.update_idletasks()

            params = self.get_simulation_params()
            if not params:
                raise RuntimeError("Invalid simulation parameters.")

            try:
                cycles = 0 if init else int(self.cycles_var.get())
            except ValueError:
                raise RuntimeError("Invalid number of cycles.")

            bounds = self.get_all_parameter_bounds()
            if not bounds:
                raise RuntimeError("Invalid parameter bounds.")
            
            simulation_type = "CCP" if self.sim_type_var.get() == 'C' else "ICP"
            self.output_text.insert(tk.END, f"Running {simulation_type} simulation in folder: {self.output_folder}\n")
            self.output_text.see(tk.END)  # Scroll to end so the message is visible
            
            
            result = self.ccp.run_simulation(
                self.output_folder.encode('utf-8'),  # NEW ARGUMENT: folder path
                opt_param.encode('utf-8'),
                self.sim_type_var.get().encode('utf-8')[0],
                cycles,
                params['target_efficiency'], params['power'], params['pressure'], params['gap_length'], params['temperature'], params['turns'], 
                bounds["power"]["min"], bounds["power"]["max"],
                bounds["pressure"]["min"], bounds["pressure"]["max"],
                bounds["gap_length"]["min"], bounds["gap_length"]["max"]
            )



            self.stop_progress()

            if result == 0:
                status = "initialisation" if init else "Simulation"
                self.output_text.insert(tk.END, f"\n{status} completed successfully\n")

                if not init:
                    self.data_dir = self.get_latest_results_dir()
                    self.output_text.insert(tk.END, f">> Found latest run folder: {self.data_dir}\n")
                    entries = self.read_heating_efficiency_log(self.data_dir)
                    self.display_efficiency_in_output(entries)
            else:
                self.output_text.insert(tk.END, f"Error code: {result}\n")

        except Exception as e:
            self.stop_progress()
            self.after(0, lambda e=e: self.output_text.insert(tk.END, f"Runtime error: {str(e)}\n"))

        finally:
            self.simulation_running = False
            self.simulation_lock.release()

    def read_heating_efficiency_log(self, run_folder):
        """
        Reads heating_efficiency_log.txt and returns list of dicts:
        If mode='N': {'cycle': int, 'efficiency': float}
        If mode='O': {'cycle': int, 'efficiency': float, 'parameter': float, 'param_name': str}
        """
        log_path = run_folder / "heating_efficiency_log.txt"
        if not log_path.exists():
            return None

        # Try to auto-detect format by checking first valid line
        entries = []
            
        with open(log_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Try full pattern first (with parameter)
                full_match = re.search(
                    r"Cycle:\s*(\d+)\s*\|\s*Efficiency:\s*([\d.]+)\s*%\s*\|\s*(\w+):\s*([\deE.+-]+)", line
                )
                if full_match:
                    cycle = int(full_match.group(1))
                    efficiency = float(full_match.group(2))
                    param_name = full_match.group(3)
                    param_value = float(full_match.group(4))
                    entries.append({
                        'cycle': cycle,
                        'efficiency': efficiency,
                        'param_name': param_name,
                        'parameter': param_value
                    })
                else:
                    # Try simple pattern (cycle and efficiency only)
                    simple_match = re.search(r"Cycle:\s*(\d+)\s*\|\s*Efficiency:\s*([\d.]+)\s*%", line)
                    if simple_match:
                        cycle = int(simple_match.group(1))
                        efficiency = float(simple_match.group(2))
                        entries.append({
                            'cycle': cycle,
                            'efficiency': efficiency
                        })

        return entries


    def display_efficiency_in_output(self, entries):
        if not entries:
            self.output_text.insert(tk.END, ">> No valid entries found in heating_efficiency_log.txt\n")
            return

        last = entries[-1]
        run_label = f">> {self.data_dir.name}:\n"
        self.output_text.insert(tk.END, run_label)
        self.output_text.insert(tk.END, f"Efficiency = {last['efficiency']:.2f}%\n")
        
        # Only display parameter info if it exists (Optimisation mode)
        if 'param_name' in last and 'parameter' in last:
            self.output_text.insert(tk.END, f"{last['param_name']} = {last['parameter']:.4g}\n")
        
        self.output_text.insert(tk.END, "\n")



if __name__ == "__main__":
    app = CCPApp()
    app.protocol("WM_DELETE_WINDOW", app.quit)
    app.mainloop()