import sys
import subprocess
from PyQt5.QtCore import Qt, pyqtSignal, QThread
from PyQt5.QtGui import QMovie, QIcon
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit, QPushButton, QCheckBox, QScrollArea, QLineEdit, QLabel, QTabWidget, QFileDialog

# GUI for ID mapper

class WebContentFetcher(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Bioinformatics Database ID Mapper')
        self.setGeometry(100, 100, 800, 600)
        icon = QIcon("dbs2.png")  # Replace with the actual path to your icon
        self.setWindowIcon(icon)

        main_layout = QVBoxLayout()

        header_label = QLabel("Bioinformatics Database ID Mapper", self)
        header_label.setAlignment(Qt.AlignCenter)
        header_label.setStyleSheet("font-size: 20px; padding: 10px; font-family: Verdana, sans-serif; color: #5A3D3D; font-weight: bold;")
        main_layout.addWidget(header_label)

        self.tab_widget = QTabWidget(self)

        # First tab
        description_tab = QWidget()
        description_layout = QVBoxLayout()
        description_text = QTextEdit()
        description_text.setPlainText("This is a description of the web content fetcher.")
        description_text.setReadOnly(True)
        description_layout.addWidget(description_text)
        description_tab.setLayout(description_layout)
        self.tab_widget.addTab(description_tab, "Description")

        # Second tab
        content_tab = QWidget()
        content_layout = QVBoxLayout()

        self.input_text = QLabel('Enter gene symbol/s or upload a file: ')
        self.input_text.setStyleSheet("font-size: 12px; font-family: Verdana, sans-serif; color: #5A3D3D; font-weight: bold;")
        self.url_input = QLineEdit()
        self.url_input.setPlaceholderText("Enter gene symbol/s seperated by comma")
        self.dblabel = QLabel('\nSelect Database/s:')
        self.dblabel.setStyleSheet("font-size: 12px; font-family: Verdana, sans-serif; color: #5A3D3D; font-weight: bold;")
        self.result_label = QLabel('\nResult:')
        self.result_label.setStyleSheet("font-size: 12px; font-family: Verdana, sans-serif; color: #5A3D3D; font-weight: bold;")
        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        self.loading_gif = QMovie("loading.gif")  # Replace "loading.gif" with the path to your GIF file
        self.loading_label = QLabel()
        self.loading_label.setMovie(self.loading_gif)

        self.fetch_button = QPushButton('Map Gene symbol/s')
        self.fetch_button.clicked.connect(self.fetch_content)
        #self.fetch_button.setFixedSize(500, 30)

        self.upload_file_button = QPushButton('Upload File')
        self.upload_file_button.clicked.connect(self.upload_file_dialog)

        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)

        self.checkbox_container = QWidget(self.scroll_area)
        self.checkbox_layout = QVBoxLayout(self.checkbox_container)
        self.checkboxes = []

        options = ["GenBank/ENA/DDBJ","UniProtKB/Swiss-Prot","PDB","GeneTree","HomoloGene","ENSEMBL","UCSC","OMIM","dbVar","SNP","CDD","Bgee","Patric","PharmGKB","eggNOG","OrthoDB","KEGG Orthology","ExpressionAtlas","TreeFam","Complex Portal","IntAct","KEGG","MINT","STRING","BRENDA","BioCyc","Reactome","BioGRID","Gene Ontology (Biological_process)","Gene Ontology (Cellular_component)","Gene Ontology (Molecular_function)"]
        for option in options:
            checkbox = QCheckBox(option)
            self.checkboxes.append(checkbox)
            self.checkbox_layout.addWidget(checkbox)

        self.scroll_area.setWidget(self.checkbox_container)
        self.scroll_area.setFixedHeight(100)

        content_layout.addWidget(self.input_text)
        content_layout.addWidget(self.url_input)
        content_layout.addWidget(self.upload_file_button)
        content_layout.addWidget(self.dblabel)
        content_layout.addWidget(self.scroll_area)
        content_layout.addWidget(self.result_label)
        content_layout.addWidget(self.result_text)
        content_layout.addWidget(self.loading_label)
        content_layout.addWidget(self.fetch_button)

        content_tab.setLayout(content_layout)
        self.tab_widget.addTab(content_tab, "Map gene symbol")

        main_layout.addWidget(self.tab_widget)
        self.setLayout(main_layout)

        # Create a thread for content fetching
        self.fetch_thread = ContentFetchThread(self)

        # Connect signals
        self.fetch_thread.contentFetched.connect(self.display_content)
        self.fetch_thread.started.connect(self.fetch_started)
        self.fetch_thread.finished.connect(self.fetch_finished)

        # Third tab
        result_tab = QWidget()
        result_layout = QVBoxLayout()
        result_text = QTextEdit()
        result_text.setPlainText("Results appear here.")
        result_text.setReadOnly(True)
        result_layout.addWidget(result_text)
        result_tab.setLayout(result_layout)
        self.tab_widget.addTab(result_tab, "Result")

    def fetch_content(self):
        input_text = self.url_input.text()
        selected_options = [checkbox.text() for checkbox in self.checkboxes if checkbox.isChecked()]

        if input_text:
            self.loading_gif.start()
            self.loading_label.show()
            self.result_text.clear()
            self.fetch_button.setEnabled(False)

            # Determine the input mode
            if input_text.startswith("http://") or input_text.startswith("https://"):
                input_mode = "URL"
            else:
                input_mode = "UploadedFile"

            self.fetch_thread.set_url_and_options(input_text, selected_options, input_mode)
            self.fetch_thread.start()
        else:
            self.result_text.setPlainText("Please enter gene symbol/s.")

    def upload_file_dialog(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Upload Text File", ",", "Text Files (*.txt);;All Files (*)", options=options)
        if file_name:
            with open(file_name, "r") as file:
                self.url_input.setText(file.read())

    def display_content(self, content):
        self.result_text.setPlainText(content)
        self.loading_gif.stop()
        self.loading_label.hide()
        self.fetch_button.setEnabled(True)

    def fetch_started(self):
        self.loading_gif.start()
        self.loading_label.show()
        self.result_text.clear()
        self.fetch_button.setEnabled(False)

    def fetch_finished(self):
        self.loading_gif.stop()
        self.loading_label.hide()
        self.fetch_button.setEnabled(True)

class ContentFetchThread(QThread):
    contentFetched = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.url = ""
        self.checkboxes = []

    def set_url_and_options(self, url, checkboxes, input_mode):
        self.url = url
        self.checkboxes = checkboxes
        self.input_mode = input_mode

    def run(self):
        content = self.get_website_content()
        self.contentFetched.emit(content)

    def get_website_content(self):
        try:
            script_path = 'fetch.py'
            arguments = [self.url] + self.checkboxes + [self.input_mode]
            arguments_to_pass = [str(arg) for arg in arguments]
            result = subprocess.check_output(['python', script_path] + arguments_to_pass, universal_newlines=True)
            return result
        except subprocess.CalledProcessError as e:
            return f"Error: {e}"



if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = WebContentFetcher()
    window.show()
    # Load and apply the external style sheet
    with open('login.qss', 'r') as qss_file:
        app.setStyleSheet(qss_file.read())
    sys.exit(app.exec_())
