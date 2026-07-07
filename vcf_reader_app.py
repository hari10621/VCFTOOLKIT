"""Windows desktop application for streaming VCF analysis."""

from __future__ import annotations

import logging
import queue
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any

from vcf_reader import VCFReader, VCFVariant, VariantFilter
from vcf_reader.filters import classify_variant

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD

    DND_AVAILABLE = True
except ImportError:  # pragma: no cover - depends on optional package
    DND_FILES = None
    TkinterDnD = None
    DND_AVAILABLE = False


APP_TITLE = "DUNIA VCF Analyzer"
PREVIEW_LIMIT = 5000


class QueueLogHandler(logging.Handler):
    """Sends log messages to the UI queue."""

    def __init__(self, events: queue.Queue[tuple[str, Any]]) -> None:
        super().__init__()
        self.events = events

    def emit(self, record: logging.LogRecord) -> None:
        self.events.put(("log", self.format(record)))


class VCFAnalyzerApp:
    """Tkinter UI for selecting, analyzing, filtering, and exporting VCF files."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1180x760")
        self.root.minsize(980, 640)

        self.events: queue.Queue[tuple[str, Any]] = queue.Queue()
        self.reader: VCFReader | None = None
        self.current_file: Path | None = None
        self.preview_variants: list[VCFVariant] = []
        self.analysis_running = False
        self.analysis_done = False
        self.export_running = False

        self.file_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Select a VCF or VCF.GZ file to begin.")
        self.progress_var = tk.DoubleVar(value=0)
        self.speed_var = tk.StringVar(value="0 variants/sec")
        self.count_var = tk.StringVar(value="0 variants")
        self.elapsed_var = tk.StringVar(value="0s")

        self._setup_logging()
        self._setup_style()
        self._build_layout()
        self._wire_drag_drop()
        self.root.after(100, self._process_events)

    def _setup_logging(self) -> None:
        handler = QueueLogHandler(self.events)
        handler.setFormatter(logging.Formatter("%(asctime)s  %(levelname)s  %(message)s", "%H:%M:%S"))
        logging.getLogger("vcf_reader").setLevel(logging.INFO)
        logging.getLogger("vcf_reader").addHandler(handler)

    def _setup_style(self) -> None:
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        self.root.configure(bg="#f6f8fb")
        style.configure("TFrame", background="#f6f8fb")
        style.configure("Surface.TFrame", background="#ffffff", borderwidth=1, relief="solid")
        style.configure("TLabel", background="#f6f8fb", foreground="#162033", font=("Segoe UI", 10))
        style.configure("Surface.TLabel", background="#ffffff", foreground="#162033", font=("Segoe UI", 10))
        style.configure("Title.TLabel", background="#f6f8fb", font=("Segoe UI Semibold", 22), foreground="#162033")
        style.configure("Subtitle.TLabel", background="#f6f8fb", font=("Segoe UI", 10), foreground="#536173")
        style.configure("Metric.TLabel", background="#ffffff", font=("Segoe UI Semibold", 18), foreground="#0f766e")
        style.configure("MetricName.TLabel", background="#ffffff", font=("Segoe UI", 9), foreground="#536173")
        style.configure("Accent.TButton", font=("Segoe UI Semibold", 10))
        style.configure("Treeview", rowheight=26, font=("Segoe UI", 9))
        style.configure("Treeview.Heading", font=("Segoe UI Semibold", 9))
        style.map("TButton", background=[("active", "#e6eef8")])

    def _build_layout(self) -> None:
        outer = ttk.Frame(self.root, padding=18)
        outer.pack(fill="both", expand=True)

        header = ttk.Frame(outer)
        header.pack(fill="x")
        ttk.Label(header, text=APP_TITLE, style="Title.TLabel").pack(side="left")
        ttk.Label(
            header,
            text="Streaming analysis for large VCF files",
            style="Subtitle.TLabel",
        ).pack(side="left", padx=(16, 0), pady=(10, 0))

        self._build_file_panel(outer)
        self._build_progress_panel(outer)
        self._build_tabs(outer)

        status = ttk.Label(outer, textvariable=self.status_var, style="Subtitle.TLabel")
        status.pack(fill="x", pady=(10, 0))

    def _build_file_panel(self, parent: ttk.Frame) -> None:
        panel = ttk.Frame(parent, style="Surface.TFrame", padding=16)
        panel.pack(fill="x", pady=(16, 12))

        left = ttk.Frame(panel, style="Surface.TFrame")
        left.pack(side="left", fill="both", expand=True)
        ttk.Label(left, text="Choose VCF file", style="Surface.TLabel", font=("Segoe UI Semibold", 12)).pack(anchor="w")
        hint = "Drop a .vcf or .vcf.gz file here, or use Browse." if DND_AVAILABLE else "Use Browse to select a .vcf or .vcf.gz file."
        ttk.Label(left, text=hint, style="Surface.TLabel").pack(anchor="w", pady=(4, 10))
        entry = ttk.Entry(left, textvariable=self.file_var)
        entry.pack(fill="x")

        actions = ttk.Frame(panel, style="Surface.TFrame")
        actions.pack(side="right", padx=(14, 0))
        ttk.Button(actions, text="Browse", command=self.choose_file).pack(fill="x", pady=(0, 8))
        self.analyze_button = ttk.Button(actions, text="Analyze", style="Accent.TButton", command=self.start_analysis)
        self.analyze_button.pack(fill="x", pady=(0, 8))
        self.cancel_button = ttk.Button(actions, text="Cancel", command=self.cancel_analysis, state="disabled")
        self.cancel_button.pack(fill="x")

        self.drop_target = panel

    def _build_progress_panel(self, parent: ttk.Frame) -> None:
        panel = ttk.Frame(parent, style="Surface.TFrame", padding=14)
        panel.pack(fill="x", pady=(0, 12))
        self.progress = ttk.Progressbar(panel, variable=self.progress_var, maximum=100)
        self.progress.pack(fill="x")

        metrics = ttk.Frame(panel, style="Surface.TFrame")
        metrics.pack(fill="x", pady=(10, 0))
        for label, variable in [
            ("Variants", self.count_var),
            ("Speed", self.speed_var),
            ("Elapsed", self.elapsed_var),
        ]:
            block = ttk.Frame(metrics, style="Surface.TFrame")
            block.pack(side="left", padx=(0, 28))
            ttk.Label(block, textvariable=variable, style="Metric.TLabel").pack(anchor="w")
            ttk.Label(block, text=label, style="MetricName.TLabel").pack(anchor="w")

    def _build_tabs(self, parent: ttk.Frame) -> None:
        self.tabs = ttk.Notebook(parent)
        self.tabs.pack(fill="both", expand=True)

        self.summary_tab = ttk.Frame(self.tabs, padding=12)
        self.header_tab = ttk.Frame(self.tabs, padding=12)
        self.stats_tab = ttk.Frame(self.tabs, padding=12)
        self.variants_tab = ttk.Frame(self.tabs, padding=12)
        self.filter_tab = ttk.Frame(self.tabs, padding=12)
        self.export_tab = ttk.Frame(self.tabs, padding=12)
        self.logs_tab = ttk.Frame(self.tabs, padding=12)

        for frame, title in [
            (self.summary_tab, "Summary"),
            (self.header_tab, "Header"),
            (self.stats_tab, "Statistics"),
            (self.variants_tab, "Variants"),
            (self.filter_tab, "Search & Filter"),
            (self.export_tab, "Export"),
            (self.logs_tab, "Logs"),
        ]:
            self.tabs.add(frame, text=title)

        self._build_summary_tab()
        self._build_header_tab()
        self._build_stats_tab()
        self._build_variants_tab()
        self._build_filter_tab()
        self._build_export_tab()
        self._build_logs_tab()

    def _build_summary_tab(self) -> None:
        self.summary_grid = ttk.Frame(self.summary_tab)
        self.summary_grid.pack(fill="x")
        self.summary_values: dict[str, tk.StringVar] = {}
        for index, name in enumerate(["File", "VCF Version", "Reference", "Samples", "Backend", "Preview Rows"]):
            variable = tk.StringVar(value="-")
            self.summary_values[name] = variable
            frame = ttk.Frame(self.summary_grid, style="Surface.TFrame", padding=12)
            frame.grid(row=index // 3, column=index % 3, sticky="ew", padx=6, pady=6)
            self.summary_grid.columnconfigure(index % 3, weight=1)
            ttk.Label(frame, text=name, style="MetricName.TLabel").pack(anchor="w")
            ttk.Label(frame, textvariable=variable, style="Metric.TLabel", wraplength=300).pack(anchor="w")

    def _build_header_tab(self) -> None:
        self.header_text = tk.Text(self.header_tab, wrap="word", height=10, font=("Consolas", 10), relief="flat")
        self.header_text.pack(fill="both", expand=True)
        self.header_text.insert("1.0", "Header details will appear after analysis.")
        self.header_text.configure(state="disabled")

    def _build_stats_tab(self) -> None:
        columns = ("metric", "value")
        self.stats_tree = ttk.Treeview(self.stats_tab, columns=columns, show="headings")
        self.stats_tree.heading("metric", text="Metric")
        self.stats_tree.heading("value", text="Value")
        self.stats_tree.column("metric", width=260, anchor="w")
        self.stats_tree.column("value", width=220, anchor="w")
        self.stats_tree.pack(fill="both", expand=True)

    def _build_variants_tab(self) -> None:
        columns = ("chromosome", "position", "id", "reference", "alternate", "qual", "filter", "gene", "type")
        self.variant_tree = ttk.Treeview(self.variants_tab, columns=columns, show="headings")
        widths = [90, 90, 130, 100, 140, 90, 110, 140, 100]
        for column, width in zip(columns, widths, strict=False):
            self.variant_tree.heading(column, text=column.title())
            self.variant_tree.column(column, width=width, anchor="w")
        self.variant_tree.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(self.variants_tab, orient="vertical", command=self.variant_tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.variant_tree.configure(yscrollcommand=scrollbar.set)

    def _build_filter_tab(self) -> None:
        controls = ttk.Frame(self.filter_tab)
        controls.pack(fill="x", pady=(0, 10))

        self.filter_chrom = tk.StringVar()
        self.filter_gene = tk.StringVar()
        self.filter_id = tk.StringVar()
        self.filter_min_qual = tk.StringVar()
        self.filter_pass_only = tk.BooleanVar()

        fields = [
            ("Chromosome", self.filter_chrom),
            ("Gene", self.filter_gene),
            ("Variant ID / rsID", self.filter_id),
            ("Min QUAL", self.filter_min_qual),
        ]
        for index, (label, variable) in enumerate(fields):
            ttk.Label(controls, text=label).grid(row=0, column=index, sticky="w", padx=4)
            ttk.Entry(controls, textvariable=variable, width=18).grid(row=1, column=index, sticky="ew", padx=4)
            controls.columnconfigure(index, weight=1)
        ttk.Checkbutton(controls, text="PASS only", variable=self.filter_pass_only).grid(row=1, column=4, padx=8)
        ttk.Button(controls, text="Run Search", command=self.run_filter_search).grid(row=1, column=5, padx=4)

        columns = ("chromosome", "position", "id", "reference", "alternate", "qual", "filter", "gene")
        self.filter_tree = ttk.Treeview(self.filter_tab, columns=columns, show="headings")
        for column in columns:
            self.filter_tree.heading(column, text=column.title())
            self.filter_tree.column(column, width=120, anchor="w")
        self.filter_tree.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(self.filter_tab, orient="vertical", command=self.filter_tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.filter_tree.configure(yscrollcommand=scrollbar.set)

    def _build_export_tab(self) -> None:
        panel = ttk.Frame(self.export_tab, style="Surface.TFrame", padding=16)
        panel.pack(anchor="nw", fill="x")
        ttk.Label(panel, text="Export filtered or full results", style="Surface.TLabel", font=("Segoe UI Semibold", 12)).pack(anchor="w")
        ttk.Label(panel, text="Choose a format and destination. Export streams from the source file.", style="Surface.TLabel").pack(anchor="w", pady=(4, 12))
        self.export_format = tk.StringVar(value="csv")
        row = ttk.Frame(panel, style="Surface.TFrame")
        row.pack(fill="x")
        ttk.Combobox(row, textvariable=self.export_format, values=["csv", "tsv", "json", "xlsx", "vcf"], width=12, state="readonly").pack(side="left")
        ttk.Button(row, text="Export All", command=lambda: self.start_export(False)).pack(side="left", padx=8)
        ttk.Button(row, text="Export Current Filter", command=lambda: self.start_export(True)).pack(side="left")

    def _build_logs_tab(self) -> None:
        self.log_text = tk.Text(self.logs_tab, wrap="word", font=("Consolas", 10), relief="flat")
        self.log_text.pack(fill="both", expand=True)
        self.log_text.configure(state="disabled")

    def _wire_drag_drop(self) -> None:
        if not DND_AVAILABLE:
            return
        self.root.drop_target_register(DND_FILES)
        self.root.dnd_bind("<<Drop>>", self._on_drop)

    def choose_file(self) -> None:
        path = filedialog.askopenfilename(
            title="Select VCF file",
            filetypes=[("VCF files", "*.vcf *.vcf.gz"), ("All files", "*.*")],
        )
        if path:
            self.file_var.set(path)
            self.status_var.set("File selected. Start analysis when ready.")

    def _on_drop(self, event: Any) -> None:
        raw = event.data.strip()
        paths = self.root.tk.splitlist(raw)
        if paths:
            self.file_var.set(paths[0])
            self.status_var.set("File dropped. Start analysis when ready.")

    def start_analysis(self) -> None:
        if self.analysis_running:
            return
        path_text = self.file_var.get().strip().strip('"')
        if not path_text:
            messagebox.showwarning(APP_TITLE, "Please select a VCF or VCF.GZ file first.")
            return
        path = Path(path_text)
        if not path.exists():
            messagebox.showerror(APP_TITLE, "The selected file does not exist.")
            return

        self.current_file = path
        self.preview_variants.clear()
        self.analysis_running = True
        self.analysis_done = False
        self.progress_var.set(0)
        self.count_var.set("0 variants")
        self.speed_var.set("0 variants/sec")
        self.elapsed_var.set("0s")
        self._clear_tree(self.variant_tree)
        self._clear_tree(self.filter_tree)
        self._clear_tree(self.stats_tree)
        self._set_text(self.header_text, "Reading header and streaming variants...")
        self.analyze_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self.status_var.set("Analyzing file...")

        thread = threading.Thread(target=self._analysis_worker, args=(path,), daemon=True)
        thread.start()
        self.root.after(250, self._tick_progress)

    def _analysis_worker(self, path: Path) -> None:
        try:
            reader = VCFReader(path)
            self.reader = reader
            reader.validate()
            header = reader.read_header()
            self.events.put(("header", header))
            for variant in reader.stream_variants(update_statistics=True):
                if len(self.preview_variants) < PREVIEW_LIMIT:
                    self.preview_variants.append(variant)
                    self.events.put(("variant", variant))
            self.events.put(("complete", reader.statistics.snapshot()))
        except Exception as exc:
            self.events.put(("error", exc))

    def cancel_analysis(self) -> None:
        if self.reader:
            self.reader.cancel()
        self.status_var.set("Cancelling analysis...")

    def _tick_progress(self) -> None:
        if self.reader and self.analysis_running:
            progress = self.reader.progress()
            self.progress_var.set(progress.percent_complete)
            self.count_var.set(f"{progress.current_variant_count:,} variants")
            self.speed_var.set(f"{progress.variants_per_second:,.0f} variants/sec")
            self.elapsed_var.set(f"{progress.elapsed_time:,.1f}s")
            self.root.after(250, self._tick_progress)

    def _process_events(self) -> None:
        try:
            while True:
                event, payload = self.events.get_nowait()
                if event == "log":
                    self._append_log(str(payload))
                elif event == "header":
                    self._show_header(payload)
                elif event == "variant":
                    self._add_variant_row(payload)
                elif event == "complete":
                    self._analysis_complete(payload)
                elif event == "error":
                    self._analysis_error(payload)
                elif event == "export_complete":
                    self._export_complete(payload)
                elif event == "export_error":
                    self._export_error(payload)
                elif event == "filter_result":
                    self.filter_tree.insert("", "end", values=self._variant_values(payload, include_type=False))
                elif event == "filter_complete":
                    self.status_var.set(f"Search complete. Showing {payload:,} matching rows.")
        except queue.Empty:
            pass
        self.root.after(100, self._process_events)

    def _show_header(self, header: Any) -> None:
        self.summary_values["File"].set(self.current_file.name if self.current_file else "-")
        self.summary_values["VCF Version"].set(header.version or "-")
        self.summary_values["Reference"].set(header.reference or "-")
        self.summary_values["Samples"].set(str(len(header.samples)))
        self.summary_values["Backend"].set(self.reader.backend.value if self.reader else "-")
        text = [
            f"VCF Version: {header.version or '-'}",
            f"Reference: {header.reference or '-'}",
            f"Samples: {', '.join(header.samples) if header.samples else '-'}",
            "",
            "INFO definitions:",
            *[f"  {key}: {value}" for key, value in header.info_definitions.items()],
            "",
            "FILTER definitions:",
            *[f"  {key}: {value}" for key, value in header.filter_definitions.items()],
            "",
            "FORMAT definitions:",
            *[f"  {key}: {value}" for key, value in header.format_definitions.items()],
            "",
            "Metadata:",
            *[f"  {key}: {value}" for key, value in header.metadata.items()],
        ]
        self._set_text(self.header_text, "\n".join(text))

    def _add_variant_row(self, variant: VCFVariant) -> None:
        self.variant_tree.insert("", "end", values=self._variant_values(variant, include_type=True))
        self.summary_values["Preview Rows"].set(f"{len(self.preview_variants):,}")

    def _analysis_complete(self, stats: Any) -> None:
        self.analysis_running = False
        self.analysis_done = True
        self.progress_var.set(100)
        self.analyze_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")
        self.status_var.set("Analysis complete.")
        self._show_stats(stats)

    def _analysis_error(self, exc: Exception) -> None:
        self.analysis_running = False
        self.analyze_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")
        self.status_var.set("Analysis failed.")
        messagebox.showerror(APP_TITLE, str(exc))
        self._append_log(f"ERROR  {exc}")

    def _show_stats(self, stats: Any) -> None:
        self._clear_tree(self.stats_tree)
        rows = [
            ("Total Variants", stats.total_variants),
            ("SNP Count", stats.snp_count),
            ("Insertion Count", stats.insertion_count),
            ("Deletion Count", stats.deletion_count),
            ("Structural Variants", stats.structural_variant_count),
            ("Transition Count", stats.transition_count),
            ("Transversion Count", stats.transversion_count),
            ("PASS Variants", stats.pass_variants),
            ("Filtered Variants", stats.filtered_variants),
            ("Average QUAL", self._format_number(stats.average_qual)),
            ("Average Depth", self._format_number(stats.average_depth)),
            ("Missing Genotypes", stats.missing_genotypes),
            ("Homozygous Count", stats.homozygous_count),
            ("Heterozygous Count", stats.heterozygous_count),
        ]
        for key, value in rows:
            self.stats_tree.insert("", "end", values=(key, value))
        for chrom, count in sorted(stats.chromosome_distribution.items()):
            self.stats_tree.insert("", "end", values=(f"Chromosome {chrom}", count))

    def run_filter_search(self) -> None:
        if not self.current_file:
            messagebox.showwarning(APP_TITLE, "Analyze a file first.")
            return
        criteria = self._current_filter()
        variant_id = self.filter_id.get().strip()
        self._clear_tree(self.filter_tree)
        self.status_var.set("Searching...")

        def worker() -> None:
            try:
                reader = VCFReader(self.current_file)
                count = 0
                for variant in reader.filter_variants(criteria):
                    if variant_id and variant.id != variant_id:
                        continue
                    self.events.put(("filter_result", variant))
                    count += 1
                    if count >= PREVIEW_LIMIT:
                        break
                self.events.put(("filter_complete", count))
                self.events.put(("log", f"Search complete. Showing {count:,} matching rows."))
            except Exception as exc:
                self.events.put(("error", exc))
        threading.Thread(target=worker, daemon=True).start()

    def start_export(self, use_filter: bool) -> None:
        if not self.current_file:
            messagebox.showwarning(APP_TITLE, "Analyze a file first.")
            return
        if self.export_running:
            return
        export_format = self.export_format.get()
        path = filedialog.asksaveasfilename(
            title="Export variants",
            defaultextension=f".{export_format}",
            filetypes=[(f"{export_format.upper()} files", f"*.{export_format}"), ("All files", "*.*")],
        )
        if not path:
            return
        self.export_running = True
        self.status_var.set("Exporting...")

        def worker() -> None:
            try:
                reader = VCFReader(self.current_file)
                reader.read_header()
                criteria = self._current_filter() if use_filter else None
                written = reader.export(path, export_format, criteria)
                self.events.put(("export_complete", (path, written)))
            except Exception as exc:
                self.events.put(("export_error", exc))

        threading.Thread(target=worker, daemon=True).start()

    def _export_complete(self, payload: tuple[str, int]) -> None:
        self.export_running = False
        path, written = payload
        self.status_var.set(f"Export complete: {written:,} variants.")
        messagebox.showinfo(APP_TITLE, f"Exported {written:,} variants to:\n{path}")

    def _export_error(self, exc: Exception) -> None:
        self.export_running = False
        self.status_var.set("Export failed.")
        messagebox.showerror(APP_TITLE, str(exc))

    def _current_filter(self) -> VariantFilter:
        min_qual_text = self.filter_min_qual.get().strip()
        min_qual = float(min_qual_text) if min_qual_text else None
        return VariantFilter(
            chromosome=self.filter_chrom.get().strip() or None,
            gene=self.filter_gene.get().strip() or None,
            min_qual=min_qual,
            pass_only=self.filter_pass_only.get(),
        )

    def _variant_values(self, variant: VCFVariant, include_type: bool) -> tuple[Any, ...]:
        gene = variant.info.get("GENE") or variant.info.get("Gene") or variant.info.get("ANN") or ""
        base = (
            variant.chromosome,
            variant.position,
            variant.id or ".",
            variant.reference,
            ",".join(variant.alternate),
            self._format_number(variant.qual),
            ";".join(variant.filter) if variant.filter else ".",
            gene,
        )
        if include_type:
            return (*base, classify_variant(variant))
        return base

    def _append_log(self, message: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def _set_text(self, widget: tk.Text, text: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        widget.configure(state="disabled")

    def _clear_tree(self, tree: ttk.Treeview) -> None:
        tree.delete(*tree.get_children())

    def _format_number(self, value: Any) -> str:
        if value is None:
            return "-"
        if isinstance(value, float):
            return f"{value:,.3f}".rstrip("0").rstrip(".")
        return str(value)


def create_root() -> tk.Tk:
    """Create a Tk root with optional drag-and-drop support."""

    if DND_AVAILABLE and TkinterDnD is not None:
        return TkinterDnD.Tk()
    return tk.Tk()


def main() -> None:
    root = create_root()
    VCFAnalyzerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
