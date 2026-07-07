"""Qt desktop application for DUNIA VCF Analyzer."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtGui import QDragEnterEvent, QDropEvent
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QStackedWidget,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextBrowser,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from vcf_reader import VCFReader, VCFVariant, VariantFilter
from vcf_reader.filters import classify_variant
from vcf_reader.inference import (
    is_clinically_significant,
    extract_clinical_metadata,
    generate_rule_based_report,
    build_gemini_prompt,
)
import json


APP_TITLE = "DUNIA VCF Analyzer"
PREVIEW_LIMIT = 5000

def get_config_path() -> Path:
    return Path.home() / ".dunia_vcf_config.json"

def load_api_key() -> str:
    try:
        path = get_config_path()
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
            return data.get("gemini_api_key", "")
    except Exception:
        pass
    return ""

def save_api_key(key: str) -> None:
    try:
        path = get_config_path()
        path.write_text(json.dumps({"gemini_api_key": key}), encoding="utf-8")
    except Exception:
        pass


class DropZone(QFrame):
    """File drop target for VCF files."""

    file_dropped = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self.setAcceptDrops(True)
        self.setObjectName("dropZone")
        layout = QVBoxLayout(self)
        title = QLabel("Drop VCF file here")
        title.setObjectName("dropTitle")
        subtitle = QLabel("Supports .vcf and .vcf.gz files. You can also use Browse.")
        subtitle.setObjectName("muted")
        title.setAlignment(Qt.AlignCenter)
        subtitle.setAlignment(Qt.AlignCenter)
        layout.addStretch()
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addStretch()

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event: QDropEvent) -> None:
        urls = event.mimeData().urls()
        if urls:
            self.file_dropped.emit(urls[0].toLocalFile())
            event.acceptProposedAction()


class AnalysisWorker(QThread):
    """Streams a VCF file in the background and collects significant variants."""

    header_ready = Signal(object)
    variant_ready = Signal(object)
    stats_ready = Signal(object)
    progress_ready = Signal(object)
    significant_ready = Signal(list)
    log_ready = Signal(str)
    failed = Signal(str)

    def __init__(self, path: Path) -> None:
        super().__init__()
        self.path = path
        self.reader: VCFReader | None = None
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True
        if self.reader:
            self.reader.cancel()

    def run(self) -> None:
        try:
            self.reader = VCFReader(self.path)
            self.log_ready.emit(f"Opened {self.path}")
            self.reader.validate()
            header = self.reader.read_header()
            self.header_ready.emit(header)
            last_progress = time.monotonic()
            
            significant_vars = []
            
            for variant in self.reader.stream_variants(update_statistics=True):
                if self._cancelled:
                    self.log_ready.emit("Analysis cancelled.")
                    return
                
                # Check for clinical significance
                if is_clinically_significant(variant) and len(significant_vars) < 2000:
                    meta = extract_clinical_metadata(variant)
                    significant_vars.append({
                        "chromosome": variant.chromosome,
                        "position": variant.position,
                        "id": variant.id or ".",
                        "ref": variant.reference,
                        "alt": ",".join(variant.alternate),
                        "qual": variant.qual,
                        "gene": meta["gene"],
                        "clnsig": meta["clnsig"],
                        "impact": meta["impact"],
                        "consequence": meta["consequence"],
                        "disease": meta["disease"],
                        "hgvsc": meta["hgvsc"],
                        "hgvsp": meta["hgvsp"]
                    })
                # Fallback: keep high quality variants as well to ensure table isn't empty
                elif len(significant_vars) < 200 and variant.qual and variant.qual >= 30:
                    meta = extract_clinical_metadata(variant)
                    significant_vars.append({
                        "chromosome": variant.chromosome,
                        "position": variant.position,
                        "id": variant.id or ".",
                        "ref": variant.reference,
                        "alt": ",".join(variant.alternate),
                        "qual": variant.qual,
                        "gene": meta["gene"],
                        "clnsig": meta["clnsig"] if meta["clnsig"] != "Not annotated" else "High Quality QC Variant",
                        "impact": meta["impact"],
                        "consequence": meta["consequence"],
                        "disease": meta["disease"],
                        "hgvsc": meta["hgvsc"],
                        "hgvsp": meta["hgvsp"]
                    })
                    
                if self.reader.progress().current_variant_count <= PREVIEW_LIMIT:
                    self.variant_ready.emit(variant)
                now = time.monotonic()
                if now - last_progress > 0.25:
                    self.progress_ready.emit(self.reader.progress())
                    last_progress = now
            self.progress_ready.emit(self.reader.progress())
            self.stats_ready.emit(self.reader.statistics.snapshot())
            self.significant_ready.emit(significant_vars)
            self.log_ready.emit("Analysis completed.")
        except Exception as exc:
            self.failed.emit(str(exc))


class FilterWorker(QThread):
    """Runs a filtered search without blocking the UI."""

    variant_ready = Signal(object)
    finished_count = Signal(int)
    failed = Signal(str)

    def __init__(self, path: Path, criteria: VariantFilter, variant_id: str | None) -> None:
        super().__init__()
        self.path = path
        self.criteria = criteria
        self.variant_id = variant_id

    def run(self) -> None:
        try:
            reader = VCFReader(self.path)
            count = 0
            for variant in reader.filter_variants(self.criteria):
                if self.variant_id and variant.id != self.variant_id:
                    continue
                self.variant_ready.emit(variant)
                count += 1
                if count >= PREVIEW_LIMIT:
                    break
            self.finished_count.emit(count)
        except Exception as exc:
            self.failed.emit(str(exc))


class ExportWorker(QThread):
    """Exports variants in the background."""

    finished_export = Signal(str, int)
    failed = Signal(str)

    def __init__(self, source: Path, output: Path, export_format: str, criteria: VariantFilter | None) -> None:
        super().__init__()
        self.source = source
        self.output = output
        self.export_format = export_format
        self.criteria = criteria

    def run(self) -> None:
        try:
            reader = VCFReader(self.source)
            reader.read_header()
            count = reader.export(self.output, self.export_format, self.criteria)
            self.finished_export.emit(str(self.output), count)
        except Exception as exc:
            self.failed.emit(str(exc))


class InferenceWorker(QThread):
    """Generates AI report in the background."""

    finished_report = Signal(str)
    failed = Signal(str)

    def __init__(self, api_key: str, model: str, prompt: str) -> None:
        super().__init__()
        self.api_key = api_key
        self.model = model
        self.prompt = prompt

    def run(self) -> None:
        try:
            from vcf_reader.inference import call_gemini_api
            report = call_gemini_api(self.api_key, self.model, self.prompt)
            self.finished_report.emit(report)
        except Exception as exc:
            self.failed.emit(str(exc))


class MainWindow(QMainWindow):
    """Main application window."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.resize(1220, 780)
        self.current_file: Path | None = None
        self.worker: AnalysisWorker | None = None
        self.filter_worker: FilterWorker | None = None
        self.export_worker: ExportWorker | None = None
        self.inference_worker: InferenceWorker | None = None
        self.preview_count = 0
        self.analysis_started = False
        self.significant_variants: list[dict[str, Any]] = []
        self.generated_report_text = ""
        self._build_ui()
        self._apply_styles()

    def _build_ui(self) -> None:
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)
        self.select_page = self._build_select_page()
        self.results_page = self._build_results_page()
        self.stack.addWidget(self.select_page)
        self.stack.addWidget(self.results_page)

    def _build_select_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(42, 36, 42, 36)
        title = QLabel(APP_TITLE)
        title.setObjectName("heroTitle")
        subtitle = QLabel("Choose a VCF file first. The app streams analysis in the background so Windows stays responsive.")
        subtitle.setObjectName("heroSubtitle")
        layout.addWidget(title)
        layout.addWidget(subtitle)

        self.drop_zone = DropZone()
        self.drop_zone.file_dropped.connect(self.set_file)
        layout.addWidget(self.drop_zone, stretch=1)

        file_row = QHBoxLayout()
        self.file_input = QLineEdit()
        self.file_input.setPlaceholderText("No file selected")
        browse = QPushButton("Browse")
        browse.clicked.connect(self.browse_file)
        self.analyze_button = QPushButton("Analyze")
        self.analyze_button.setObjectName("primaryButton")
        self.analyze_button.clicked.connect(self.start_analysis)
        file_row.addWidget(self.file_input, stretch=1)
        file_row.addWidget(browse)
        file_row.addWidget(self.analyze_button)
        layout.addLayout(file_row)
        return page

    def _build_results_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(18, 18, 18, 18)

        top = QHBoxLayout()
        self.file_label = QLabel("No file")
        self.file_label.setObjectName("sectionTitle")
        self.status_label = QLabel("Waiting for file")
        self.status_label.setObjectName("statusLabel")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_analysis)
        back = QPushButton("Choose Another File")
        back.clicked.connect(lambda: self.stack.setCurrentWidget(self.select_page))
        top.addWidget(self.file_label, stretch=1)
        top.addWidget(self.status_label)
        top.addWidget(self.cancel_button)
        top.addWidget(back)
        layout.addLayout(top)

        progress_row = QGridLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.variant_metric = QLabel("0 variants")
        self.speed_metric = QLabel("0 variants/sec")
        self.elapsed_metric = QLabel("0s")
        progress_row.addWidget(self.progress_bar, 0, 0, 1, 3)
        progress_row.addWidget(self.variant_metric, 1, 0)
        progress_row.addWidget(self.speed_metric, 1, 1)
        progress_row.addWidget(self.elapsed_metric, 1, 2)
        layout.addLayout(progress_row)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, stretch=1)
        self._build_summary_tab()
        self._build_header_tab()
        self._build_stats_tab()
        self._build_variants_tab()
        self._build_inference_tab()
        self._build_filter_tab()
        self._build_export_tab()
        self._build_logs_tab()
        return page

    def _build_summary_tab(self) -> None:
        tab = QWidget()
        layout = QGridLayout(tab)
        self.summary_labels: dict[str, QLabel] = {}
        for index, name in enumerate(["File", "VCF Version", "Reference", "Samples", "Backend", "Preview Rows"]):
            frame = QFrame()
            frame.setObjectName("metricBox")
            box = QVBoxLayout(frame)
            label = QLabel(name)
            label.setObjectName("muted")
            value = QLabel("-")
            value.setObjectName("metricValue")
            value.setWordWrap(True)
            box.addWidget(label)
            box.addWidget(value)
            layout.addWidget(frame, index // 3, index % 3)
            self.summary_labels[name] = value
        self.tabs.addTab(tab, "Summary")

    def _build_header_tab(self) -> None:
        self.header_text = QTextEdit()
        self.header_text.setReadOnly(True)
        self.header_text.setPlainText("Header details will appear after analysis.")
        self.tabs.addTab(self.header_text, "Header")

    def _build_stats_tab(self) -> None:
        self.stats_table = self._new_table(["Metric", "Value"])
        self.tabs.addTab(self.stats_table, "Statistics")

    def _build_variants_tab(self) -> None:
        self.variant_table = self._new_table(["Chromosome", "Position", "ID", "Reference", "Alternate", "QUAL", "Filter", "Gene", "Type"])
        self.tabs.addTab(self.variant_table, "Variants")

    def _build_inference_tab(self) -> None:
        tab = QWidget()
        main_layout = QHBoxLayout(tab)
        
        # Splitter to separate Left Panel (Table/Details) and Right Panel (AI Report)
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # LEFT PANEL
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        
        # Dashboard cards grid
        cards_widget = QWidget()
        cards_layout = QGridLayout(cards_widget)
        cards_layout.setContentsMargins(0, 0, 0, 0)
        
        self.kpi_labels = {}
        kpi_names = [
            ("Pathogenic", "pathogenicCard", "0", "#b91c1c"),
            ("Likely Pathogenic", "likelyCard", "0", "#c2410c"),
            ("High Impact Consequence", "impactCard", "0", "#6b21a8"),
            ("VUS / Uncertain", "vusCard", "0", "#1e3a8a"),
            ("Unique Affected Genes", "genesCard", "0", "#115e59")
        ]
        
        for index, (name, obj_id, val, color) in enumerate(kpi_names):
            frame = QFrame()
            frame.setObjectName("kpiCard")
            
            # Direct style setting for clean cross-platform aesthetics
            bg_colors = {
                "pathogenicCard": ("#fef2f2", "#ef4444"),
                "likelyCard": ("#fff7ed", "#f97316"),
                "impactCard": ("#faf5ff", "#a855f7"),
                "vusCard": ("#eff6ff", "#3b82f6"),
                "genesCard": ("#f0fdfa", "#14b8a6")
            }
            bg, border = bg_colors[obj_id]
            frame.setStyleSheet(f"background: {bg}; border: 1px solid #d8dee9; border-left: 5px solid {border}; border-radius: 8px; padding: 6px;")
            
            box = QVBoxLayout(frame)
            box.setContentsMargins(8, 8, 8, 8)
            
            lbl_title = QLabel(name)
            lbl_title.setObjectName("kpiTitle")
            lbl_title.setStyleSheet("font-size: 8pt; color: #5c6678; font-weight: 600; background: transparent;")
            
            lbl_val = QLabel(val)
            lbl_val.setObjectName("kpiValue")
            lbl_val.setStyleSheet(f"font-size: 18px; font-weight: 700; color: {color}; background: transparent;")
            
            box.addWidget(lbl_title)
            box.addWidget(lbl_val)
            
            cards_layout.addWidget(frame, index // 3, index % 3)
            self.kpi_labels[name] = lbl_val
            
        left_layout.addWidget(cards_widget)
        
        # Significant Variants Table
        table_lbl = QLabel("Identified Significant / High-Impact Variants")
        table_lbl.setStyleSheet("font-weight: bold; font-size: 11pt; margin-top: 10px; color: #1e293b;")
        left_layout.addWidget(table_lbl)
        
        self.inference_table = self._new_table(["Gene", "Chr", "Pos", "rsID", "Consequence", "ClinVar Significance", "Impact"])
        self.inference_table.itemSelectionChanged.connect(self.on_inference_variant_selected)
        left_layout.addWidget(self.inference_table, stretch=2)
        
        # Variant Details Box
        self.detail_group = QGroupBox("Selected Variant Details")
        self.detail_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #cbd6e2;
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 12px;
                background-color: white;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 3px 0 3px;
                color: #0f766e;
            }
        """)
        detail_layout = QGridLayout(self.detail_group)
        detail_layout.setContentsMargins(12, 12, 12, 12)
        
        self.detail_fields = {}
        detail_labels = [
            ("Gene / Transcript:", 0, 0),
            ("Chr / Pos:", 0, 2),
            ("HGVS.c (Coding):", 1, 0),
            ("HGVS.p (Protein):", 1, 2),
            ("Consequence / Impact:", 2, 0),
            ("ClinVar Significance:", 2, 2),
            ("Associated Disease:", 3, 0)
        ]
        
        for name, row, col in detail_labels:
            lbl = QLabel(name)
            lbl.setStyleSheet("font-weight: bold; color: #475569; background: transparent;")
            val = QLabel("-")
            val.setWordWrap(True)
            val.setStyleSheet("color: #0f172a; background: transparent;")
            detail_layout.addWidget(lbl, row, col)
            if name == "Associated Disease:":
                detail_layout.addWidget(val, row, col + 1, 1, 3)
            else:
                detail_layout.addWidget(val, row, col + 1)
            self.detail_fields[name] = val
            
        left_layout.addWidget(self.detail_group)
        
        # RIGHT PANEL (AI Report)
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        ai_box = QGroupBox("AI Clinical Interpretation")
        ai_box.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #cbd6e2;
                border-radius: 8px;
                padding-top: 12px;
                background-color: white;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 3px 0 3px;
                color: #0f766e;
            }
        """)
        ai_layout = QVBoxLayout(ai_box)
        ai_layout.setContentsMargins(12, 12, 12, 12)
        
        # Input row for API Key
        key_layout = QHBoxLayout()
        key_lbl = QLabel("Gemini API Key:")
        key_lbl.setStyleSheet("font-weight: bold; color: #475569; background: transparent;")
        self.api_key_input = QLineEdit()
        self.api_key_input.setEchoMode(QLineEdit.Password)
        self.api_key_input.setPlaceholderText("Enter Gemini API Key (Optional)...")
        # Load saved API key
        self.api_key_input.setText(load_api_key())
        
        self.save_key_check = QCheckBox("Save Key")
        self.save_key_check.setChecked(bool(load_api_key()))
        self.save_key_check.setStyleSheet("background: transparent;")
        
        key_layout.addWidget(key_lbl)
        key_layout.addWidget(self.api_key_input, stretch=1)
        key_layout.addWidget(self.save_key_check)
        ai_layout.addLayout(key_layout)
        
        # Option row
        opt_layout = QHBoxLayout()
        model_lbl = QLabel("Model:")
        model_lbl.setStyleSheet("font-weight: bold; color: #475569; background: transparent;")
        self.ai_model_select = QComboBox()
        self.ai_model_select.addItems(["gemini-2.5-flash", "gemini-2.5-pro"])
        
        type_lbl = QLabel("Report:")
        type_lbl.setStyleSheet("font-weight: bold; color: #475569; background: transparent;")
        self.report_type_select = QComboBox()
        self.report_type_select.addItems([
            "Clinical Executive Summary", 
            "Detailed Pathogenicity Analysis", 
            "Therapeutic and Clinical Trial Actionability"
        ])
        
        opt_layout.addWidget(model_lbl)
        opt_layout.addWidget(self.ai_model_select)
        opt_layout.addWidget(type_lbl)
        opt_layout.addWidget(self.report_type_select, stretch=1)
        ai_layout.addLayout(opt_layout)
        
        # Action buttons
        btn_layout = QHBoxLayout()
        self.btn_gen_ai = QPushButton("Generate AI Report")
        self.btn_gen_ai.setObjectName("primaryButton")
        self.btn_gen_ai.setStyleSheet("""
            QPushButton#primaryButton {
                background: #0f766e;
                color: white;
                border: 1px solid #0f766e;
                font-weight: bold;
                border-radius: 6px;
                padding: 8px 14px;
            }
            QPushButton#primaryButton:hover {
                background: #0d5e58;
            }
        """)
        self.btn_gen_ai.clicked.connect(self.run_ai_inference)
        
        self.btn_gen_offline = QPushButton("Generate Offline Report")
        self.btn_gen_offline.clicked.connect(self.run_offline_inference)
        
        btn_layout.addWidget(self.btn_gen_ai)
        btn_layout.addWidget(self.btn_gen_offline)
        ai_layout.addLayout(btn_layout)
        
        # Report text area
        self.report_viewer = QTextBrowser()
        self.report_viewer.setOpenExternalLinks(True)
        self.report_viewer.setPlaceholderText("The clinical report will appear here. If no Gemini API Key is provided, you can generate an Offline report using annotations extracted from the VCF file.")
        self.report_viewer.setStyleSheet("border: 1px solid #cbd6e2; border-radius: 6px; padding: 10px; background-color: #f8fafc; color: #0f172a;")
        ai_layout.addWidget(self.report_viewer, stretch=1)
        
        # Export button
        export_layout = QHBoxLayout()
        self.btn_export_report = QPushButton("Export Report")
        self.btn_export_report.clicked.connect(self.export_clinical_report)
        self.btn_export_report.setEnabled(False)
        export_layout.addStretch()
        export_layout.addWidget(self.btn_export_report)
        ai_layout.addLayout(export_layout)
        
        right_layout.addWidget(ai_box)
        
        # Add to splitter
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        
        self.tabs.addTab(tab, "Clinical Inference")

    def _build_filter_tab(self) -> None:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        form_row = QHBoxLayout()
        self.chrom_filter = QLineEdit()
        self.gene_filter = QLineEdit()
        self.id_filter = QLineEdit()
        self.qual_filter = QLineEdit()
        self.pass_filter = QCheckBox("PASS only")
        for widget, placeholder in [
            (self.chrom_filter, "Chromosome"),
            (self.gene_filter, "Gene"),
            (self.id_filter, "Variant ID / rsID"),
            (self.qual_filter, "Min QUAL"),
        ]:
            widget.setPlaceholderText(placeholder)
            form_row.addWidget(widget)
        form_row.addWidget(self.pass_filter)
        search = QPushButton("Run Search")
        search.clicked.connect(self.run_filter_search)
        form_row.addWidget(search)
        layout.addLayout(form_row)
        self.filter_table = self._new_table(["Chromosome", "Position", "ID", "Reference", "Alternate", "QUAL", "Filter", "Gene"])
        layout.addWidget(self.filter_table)
        self.tabs.addTab(tab, "Search & Filter")

    def _build_export_tab(self) -> None:
        tab = QWidget()
        layout = QFormLayout(tab)
        self.export_format = QComboBox()
        self.export_format.addItems(["csv", "tsv", "json", "xlsx", "vcf"])
        export_all = QPushButton("Export All")
        export_all.clicked.connect(lambda: self.start_export(False))
        export_filtered = QPushButton("Export Current Filter")
        export_filtered.clicked.connect(lambda: self.start_export(True))
        layout.addRow("Format", self.export_format)
        layout.addRow(export_all)
        layout.addRow(export_filtered)
        self.tabs.addTab(tab, "Export")

    def _build_logs_tab(self) -> None:
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.tabs.addTab(self.log_text, "Logs")

    def _new_table(self, columns: list[str]) -> QTableWidget:
        table = QTableWidget(0, len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        return table

    def browse_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select VCF file", "", "VCF files (*.vcf *.vcf.gz);;All files (*)")
        if path:
            self.set_file(path, auto_start=True)

    def set_file(self, path: str, auto_start: bool = True) -> None:
        self.file_input.setText(path)
        self.current_file = Path(path)
        if auto_start:
            self.start_analysis()

    def start_analysis(self) -> None:
        if self.worker and self.worker.isRunning():
            return
        path = Path(self.file_input.text().strip())
        if not path.exists():
            QMessageBox.warning(self, APP_TITLE, "Please select a valid VCF or VCF.GZ file first.")
            return
        self.current_file = path
        self.analysis_started = True
        self.stack.setCurrentWidget(self.results_page)
        self.file_label.setText(path.name)
        self.status_label.setText("Analyzing...")
        self.preview_count = 0
        self._clear_results()
        self.tabs.setCurrentIndex(0)
        self.cancel_button.setEnabled(True)
        self.worker = AnalysisWorker(path)
        self.worker.header_ready.connect(self.show_header)
        self.worker.variant_ready.connect(self.add_variant)
        self.worker.stats_ready.connect(self.show_stats)
        self.worker.progress_ready.connect(self.show_progress)
        self.worker.significant_ready.connect(self.show_inference)
        self.worker.log_ready.connect(self.add_log)
        self.worker.failed.connect(self.analysis_failed)
        self.worker.finished.connect(lambda: self.cancel_button.setEnabled(False))
        self.worker.start()

    def cancel_analysis(self) -> None:
        if self.worker:
            self.worker.cancel()
            self.cancel_button.setEnabled(False)
            self.status_label.setText("Cancelling...")

    def show_inference(self, significant_vars: list[dict[str, Any]]) -> None:
        self.significant_variants = significant_vars
        
        # Calculate statistics
        pathogenic = sum(1 for v in significant_vars if "pathogenic" in v["clnsig"].lower() and "likely" not in v["clnsig"].lower())
        likely = sum(1 for v in significant_vars if "likely pathogenic" in v["clnsig"].lower())
        high_impact = sum(1 for v in significant_vars if v["impact"].upper() == "HIGH")
        vus = sum(1 for v in significant_vars if "uncertain" in v["clnsig"].lower() or "vus" in v["clnsig"].lower())
        unique_genes = len(set(v["gene"] for v in significant_vars if v["gene"] != "Unknown"))
        
        self.kpi_labels["Pathogenic"].setText(str(pathogenic))
        self.kpi_labels["Likely Pathogenic"].setText(str(likely))
        self.kpi_labels["High Impact Consequence"].setText(str(high_impact))
        self.kpi_labels["VUS / Uncertain"].setText(str(vus))
        self.kpi_labels["Unique Affected Genes"].setText(str(unique_genes))
        
        # Populate table
        self.inference_table.setRowCount(0)
        for v in significant_vars:
            row = self.inference_table.rowCount()
            self.inference_table.insertRow(row)
            
            # ["Gene", "Chr", "Pos", "rsID", "Consequence", "ClinVar Significance", "Impact"]
            self.inference_table.setItem(row, 0, QTableWidgetItem(str(v["gene"])))
            self.inference_table.setItem(row, 1, QTableWidgetItem(str(v["chromosome"])))
            self.inference_table.setItem(row, 2, QTableWidgetItem(str(v["position"])))
            self.inference_table.setItem(row, 3, QTableWidgetItem(str(v["id"])))
            self.inference_table.setItem(row, 4, QTableWidgetItem(str(v["consequence"])))
            self.inference_table.setItem(row, 5, QTableWidgetItem(str(v["clnsig"])))
            
            # Format impact column nicely
            impact_item = QTableWidgetItem(str(v["impact"]))
            if v["impact"].upper() == "HIGH":
                impact_item.setForeground(Qt.red)
            elif v["impact"].upper() == "MODERATE":
                impact_item.setForeground(Qt.darkYellow)
            self.inference_table.setItem(row, 6, impact_item)
            
        if significant_vars:
            self.inference_table.selectRow(0)
            self.on_inference_variant_selected()

    def on_inference_variant_selected(self) -> None:
        selected = self.inference_table.currentRow()
        if selected < 0 or selected >= len(self.significant_variants):
            return
            
        v = self.significant_variants[selected]
        
        self.detail_fields["Gene / Transcript:"].setText(f"{v['gene']}   (Transcript: {v['hgvsp'] if v['hgvsp'] else 'N/A'})")
        self.detail_fields["Chr / Pos:"].setText(f"Chromosome {v['chromosome']} : {v['position']:,}")
        self.detail_fields["HGVS.c (Coding):"].setText(v["hgvsc"] or "-")
        self.detail_fields["HGVS.p (Protein):"].setText(v["hgvsp"] or "N/A")
        self.detail_fields["Consequence / Impact:"].setText(f"{v['consequence']}   [Impact: {v['impact']}]")
        self.detail_fields["ClinVar Significance:"].setText(v["clnsig"])
        self.detail_fields["Associated Disease:"].setText(v["disease"])

    def run_ai_inference(self) -> None:
        if not self.significant_variants:
            QMessageBox.warning(self, APP_TITLE, "No significant variants loaded. Please analyze a VCF file first.")
            return
            
        api_key = self.api_key_input.text().strip()
        if not api_key:
            QMessageBox.warning(self, APP_TITLE, "Please enter your Gemini API Key first.\nOr use 'Generate Offline Report' to get an instant analysis without a key.")
            return
            
        # Manage saving API Key
        if self.save_key_check.isChecked():
            save_api_key(api_key)
        else:
            # Delete config file if unchecked
            try:
                get_config_path().unlink(missing_ok=True)
            except Exception:
                pass
                
        self.btn_gen_ai.setEnabled(False)
        self.btn_gen_offline.setEnabled(False)
        self.report_viewer.setPlainText("Contacting Google Gemini API and generating clinical report...\nThis may take up to 20 seconds. Please wait.")
        
        model = self.ai_model_select.currentText()
        report_type = self.report_type_select.currentText()
        prompt = build_gemini_prompt(self.significant_variants, self.current_file.name if self.current_file else "unknown.vcf", report_type)
        
        self.inference_worker = InferenceWorker(api_key, model, prompt)
        self.inference_worker.finished_report.connect(self.on_ai_report_finished)
        self.inference_worker.failed.connect(self.on_ai_report_failed)
        self.inference_worker.start()
        
    def on_ai_report_finished(self, report: str) -> None:
        self.generated_report_text = report
        self.report_viewer.setMarkdown(report)
        self.btn_export_report.setEnabled(True)
        self.btn_gen_ai.setEnabled(True)
        self.btn_gen_offline.setEnabled(True)
        self.add_log("AI clinical interpretation report generated.")
        
    def on_ai_report_failed(self, error: str) -> None:
        self.report_viewer.setPlainText(f"Failed to generate AI report.\n\nError details:\n{error}")
        QMessageBox.critical(self, APP_TITLE, f"AI report generation failed:\n{error}")
        self.btn_gen_ai.setEnabled(True)
        self.btn_gen_offline.setEnabled(True)
        
    def run_offline_inference(self) -> None:
        if not self.significant_variants:
            QMessageBox.warning(self, APP_TITLE, "No variants loaded. Please analyze a VCF file first.")
            return
            
        report = generate_rule_based_report(self.significant_variants, self.current_file.name if self.current_file else "unknown.vcf")
        self.generated_report_text = report
        self.report_viewer.setMarkdown(report)
        self.btn_export_report.setEnabled(True)
        self.add_log("Offline clinical report generated.")
        
    def export_clinical_report(self) -> None:
        if not self.generated_report_text:
            return
            
        path, selected_filter = QFileDialog.getSaveFileName(
            self, 
            "Export Clinical Report", 
            "Genomic_Clinical_Report.md", 
            "Markdown files (*.md);;HTML files (*.html);;All files (*)"
        )
        if not path:
            return
            
        try:
            out_path = Path(path)
            if out_path.suffix.lower() == ".html":
                out_path.write_text(self.report_viewer.toHtml(), encoding="utf-8")
            else:
                out_path.write_text(self.generated_report_text, encoding="utf-8")
            QMessageBox.information(self, APP_TITLE, f"Report exported successfully to:\n{path}")
            self.add_log(f"Report exported to {path}")
        except Exception as exc:
            QMessageBox.critical(self, APP_TITLE, f"Failed to export report: {exc}")

    def show_header(self, header: Any) -> None:
        self.status_label.setText("Header parsed")
        self.summary_labels["File"].setText(self.current_file.name if self.current_file else "-")
        self.summary_labels["VCF Version"].setText(header.version or "-")
        self.summary_labels["Reference"].setText(header.reference or "-")
        self.summary_labels["Samples"].setText(str(len(header.samples)))
        self.summary_labels["Backend"].setText(self.worker.reader.backend.value if self.worker and self.worker.reader else "-")
        lines = [
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
        ]
        self.header_text.setPlainText("\n".join(lines))

    def add_variant(self, variant: VCFVariant) -> None:
        self.preview_count += 1
        self.summary_labels["Preview Rows"].setText(f"{self.preview_count:,}")
        self._append_row(self.variant_table, self._variant_values(variant, include_type=True))
        if self.preview_count == 1:
            self.status_label.setText("Variants found")

    def show_progress(self, progress: Any) -> None:
        self.progress_bar.setValue(int(progress.percent_complete))
        self.variant_metric.setText(f"{progress.current_variant_count:,} variants")
        self.speed_metric.setText(f"{progress.variants_per_second:,.0f} variants/sec")
        self.elapsed_metric.setText(f"{progress.elapsed_time:,.1f}s")

    def show_stats(self, stats: Any) -> None:
        self.progress_bar.setValue(100)
        self.status_label.setText(f"Complete - {stats.total_variants:,} variants analyzed")
        self.stats_table.setRowCount(0)
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
            ("Average QUAL", self._format_value(stats.average_qual)),
            ("Average Depth", self._format_value(stats.average_depth)),
            ("Missing Genotypes", stats.missing_genotypes),
            ("Homozygous Count", stats.homozygous_count),
            ("Heterozygous Count", stats.heterozygous_count),
        ]
        for row in rows:
            self._append_row(self.stats_table, row)
        for chrom, count in sorted(stats.chromosome_distribution.items()):
            self._append_row(self.stats_table, (f"Chromosome {chrom}", count))
        self.tabs.setCurrentWidget(self.stats_table)
        if stats.total_variants == 0:
            self.add_log("Analysis complete, but no variant rows were found after the header.")
            QMessageBox.information(self, APP_TITLE, "Analysis finished, but no variant rows were found.")
        else:
            self.add_log(f"Analysis complete. {stats.total_variants:,} variants analyzed. Preview rows: {self.preview_count:,}.")

    def analysis_failed(self, message: str) -> None:
        self.cancel_button.setEnabled(False)
        self.status_label.setText("Failed")
        self.add_log(f"ERROR: {message}")
        QMessageBox.critical(self, APP_TITLE, message)

    def run_filter_search(self) -> None:
        if not self.current_file:
            QMessageBox.warning(self, APP_TITLE, "Analyze a file first.")
            return
        self.filter_table.setRowCount(0)
        self.filter_worker = FilterWorker(self.current_file, self._current_filter(), self.id_filter.text().strip() or None)
        self.filter_worker.variant_ready.connect(lambda variant: self._append_row(self.filter_table, self._variant_values(variant, include_type=False)))
        self.filter_worker.finished_count.connect(lambda count: self.add_log(f"Search complete. Showing {count:,} rows."))
        self.filter_worker.failed.connect(self.analysis_failed)
        self.filter_worker.start()

    def start_export(self, use_filter: bool) -> None:
        if not self.current_file:
            QMessageBox.warning(self, APP_TITLE, "Analyze a file first.")
            return
        fmt = self.export_format.currentText()
        path, _ = QFileDialog.getSaveFileName(self, "Export variants", f"variants.{fmt}", f"{fmt.upper()} files (*.{fmt});;All files (*)")
        if not path:
            return
        criteria = self._current_filter() if use_filter else None
        self.export_worker = ExportWorker(self.current_file, Path(path), fmt, criteria)
        self.export_worker.finished_export.connect(self.export_finished)
        self.export_worker.failed.connect(self.analysis_failed)
        self.export_worker.start()
        self.add_log("Export started.")

    def export_finished(self, path: str, count: int) -> None:
        self.add_log(f"Export complete: {count:,} variants -> {path}")
        QMessageBox.information(self, APP_TITLE, f"Exported {count:,} variants to:\n{path}")

    def _current_filter(self) -> VariantFilter:
        qual = self.qual_filter.text().strip()
        try:
            min_qual = float(qual) if qual else None
        except ValueError:
            QMessageBox.warning(self, APP_TITLE, "Min QUAL must be a number.")
            min_qual = None
        return VariantFilter(
            chromosome=self.chrom_filter.text().strip() or None,
            gene=self.gene_filter.text().strip() or None,
            min_qual=min_qual,
            pass_only=self.pass_filter.isChecked(),
        )

    def _variant_values(self, variant: VCFVariant, include_type: bool) -> tuple[Any, ...]:
        gene = variant.info.get("GENE") or variant.info.get("Gene") or variant.info.get("ANN") or ""
        values = (
            variant.chromosome,
            variant.position,
            variant.id or ".",
            variant.reference,
            ",".join(variant.alternate),
            self._format_value(variant.qual),
            ";".join(variant.filter) if variant.filter else ".",
            gene,
        )
        if include_type:
            return (*values, classify_variant(variant))
        return values

    def _append_row(self, table: QTableWidget, values: tuple[Any, ...]) -> None:
        row = table.rowCount()
        table.insertRow(row)
        for column, value in enumerate(values):
            table.setItem(row, column, QTableWidgetItem(str(value)))

    def _clear_results(self) -> None:
        self.progress_bar.setValue(0)
        self.variant_metric.setText("0 variants")
        self.speed_metric.setText("0 variants/sec")
        self.elapsed_metric.setText("0s")
        self.variant_table.setRowCount(0)
        self.filter_table.setRowCount(0)
        self.stats_table.setRowCount(0)
        # Clear inference elements
        self.significant_variants = []
        if hasattr(self, "inference_table"):
            self.inference_table.setRowCount(0)
        if hasattr(self, "report_viewer"):
            self.report_viewer.clear()
        self.generated_report_text = ""
        if hasattr(self, "btn_export_report"):
            self.btn_export_report.setEnabled(False)
        if hasattr(self, "kpi_labels"):
            for label in self.kpi_labels.values():
                label.setText("0")
        if hasattr(self, "detail_fields"):
            for label in self.detail_fields.values():
                label.setText("-")
        
        self.log_text.clear()
        self.header_text.setPlainText("Reading header and streaming variants...")
        self.add_log("Analysis started.")
        for label in self.summary_labels.values():
            label.setText("-")

    def add_log(self, message: str) -> None:
        self.log_text.append(message)

    def _format_value(self, value: Any) -> str:
        if value is None:
            return "-"
        if isinstance(value, float):
            return f"{value:,.3f}".rstrip("0").rstrip(".")
        return str(value)

    def _apply_styles(self) -> None:
        self.setStyleSheet(
            """
            QWidget {
                background: #f5f7fb;
                color: #172033;
                font-family: Segoe UI;
                font-size: 10pt;
            }
            QFrame#dropZone {
                background: #ffffff;
                border: 2px dashed #5a7fb7;
                border-radius: 8px;
                min-height: 260px;
            }
            QLabel#heroTitle {
                font-size: 28px;
                font-weight: 700;
            }
            QLabel#heroSubtitle, QLabel#muted {
                color: #5c6678;
            }
            QLabel#dropTitle {
                font-size: 22px;
                font-weight: 700;
                color: #1f5f99;
                background: #ffffff;
            }
            QLabel#sectionTitle {
                font-size: 18px;
                font-weight: 700;
            }
            QLabel#statusLabel {
                color: #0f766e;
                font-weight: 700;
                padding: 6px 10px;
                background: #e7f6f3;
                border: 1px solid #b8ded7;
                border-radius: 6px;
            }
            QLabel#metricValue {
                font-size: 20px;
                font-weight: 700;
                color: #0f766e;
            }
            QFrame#metricBox {
                background: #ffffff;
                border: 1px solid #d8dee9;
                border-radius: 8px;
                padding: 10px;
            }
            QPushButton {
                background: #ffffff;
                border: 1px solid #b9c5d6;
                border-radius: 6px;
                padding: 8px 14px;
            }
            QPushButton:hover {
                background: #eef4fb;
            }
            QPushButton#primaryButton {
                background: #0f766e;
                color: white;
                border: 1px solid #0f766e;
                font-weight: 700;
            }
            QLineEdit, QTextEdit, QComboBox {
                background: white;
                border: 1px solid #c8d2df;
                border-radius: 6px;
                padding: 7px;
            }
            QTableWidget {
                background: white;
                alternate-background-color: #f0f5f9;
                gridline-color: #d9e1ea;
            }
            QHeaderView::section {
                background: #e8eef6;
                padding: 6px;
                border: 1px solid #cbd6e2;
                font-weight: 700;
            }
            """
        )


def main() -> int:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
