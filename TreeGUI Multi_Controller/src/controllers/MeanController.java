package controllers;

import java.awt.Desktop;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.ResourceBundle;

import centroid.CentroidMain;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.scene.control.Button;
import javafx.scene.control.ListView;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleButton;
import javafx.scene.control.ToggleGroup;
import javafx.scene.input.DragEvent;
import javafx.scene.input.Dragboard;
import javafx.scene.input.TransferMode;
import javafx.scene.text.Text;

public class MeanController {
	String options;
	String inFile;
	String outFile;
	String command;
	File folder;
	ArrayList<File> meanfiles = new ArrayList<File>();
	ObservableList<String> mean_filenameObv = FXCollections
			.observableArrayList("Input Files Here");
	int iterations = 0;
	File log;

	@FXML
	// ResourceBundle that was given to the FXMLLoader
	private ResourceBundle resources;

	@FXML
	// URL location of the FXML file that was given to the FXMLLoader
	private URL location;

	@FXML
	// fx:id="eValue_field"
	private TextField eValue_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="rand_radio"
	private RadioButton rand_radio; // Value injected by FXMLLoader

	@FXML
	// fx:id="randperm_radio"
	private RadioButton randperm_radio; // Value injected by FXMLLoader

	@FXML
	// fx:id="iterations_min"
	private Text iterations_min; // Value injected by FXMLLoader

	@FXML
	// fx:id="iterations_slide"
	private Slider iterations_slide; // Value injected by FXMLLoader

	@FXML
	// fx:id="perm_label"
	private Text perm_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanOut_label"
	private Text meanOut_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="perm_toggle"
	private ToggleButton perm_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="iterations_label"
	private Text iterations_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="eValue_label"
	private Text eValue_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="interval_label"
	private Text interval_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="cauchy_label"
	private Text cauchy_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="interval_field"
	private TextField interval_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="iterations_max"
	private Text iterations_max; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanOut_field"
	private TextField meanOut_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="cauchy_field"
	private TextField cauchy_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="eFactor_field"
	private TextField eFactor_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="algLabel"
	private Text algLabel; // Value injected by FXMLLoader

	@FXML
	// fx:id="eFactorLabel"
	private Text eFactorLabel; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanRooted_toggle"
	private ToggleButton meanRooted_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="sturmSubmit_button"
	private Button sturmSubmit_button; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanRooted_label"
	private Text meanRooted_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanSummary_label"
	private Text meanSummary_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="mean_file_list"
	private ListView<String> mean_file_list; // Value injected by FXMLLoader

	@FXML
	// fx:id="meanTabPane"
	private SplitPane meanTabPane; // Value injected by FXMLLoader

	@FXML
	// This method is called by the FXMLLoader when initialization is complete
	void initialize() {
		assert eValue_field != null : "fx:id=\"eValue_field\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert rand_radio != null : "fx:id=\"rand_radio\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert randperm_radio != null : "fx:id=\"randperm_radio\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert iterations_min != null : "fx:id=\"iterations_min\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert iterations_slide != null : "fx:id=\"iterations_slide\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert perm_label != null : "fx:id=\"perm_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanOut_label != null : "fx:id=\"meanOut_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert perm_toggle != null : "fx:id=\"perm_toggle\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert iterations_label != null : "fx:id=\"iterations_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert eValue_label != null : "fx:id=\"eValue_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert interval_label != null : "fx:id=\"interval_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert cauchy_label != null : "fx:id=\"cauchy_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert interval_field != null : "fx:id=\"interval_field\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert iterations_max != null : "fx:id=\"iterations_max\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanOut_field != null : "fx:id=\"meanOut_field\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert cauchy_field != null : "fx:id=\"cauchy_field\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert eFactor_field != null : "fx:id=\"eFactor_field\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert algLabel != null : "fx:id=\"algLabel\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert eFactorLabel != null : "fx:id=\"eFactorLabel\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanRooted_toggle != null : "fx:id=\"meanRooted_toggle\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert sturmSubmit_button != null : "fx:id=\"sturmSubmit_button\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanRooted_label != null : "fx:id=\"meanRooted_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanSummary_label != null : "fx:id=\"meanSummary_label\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert mean_file_list != null : "fx:id=\"mean_file_list\" was not injected: check your FXML file 'MeanOptions.fxml'.";
		assert meanTabPane != null : "fx:id=\"meanTabPane\" was not injected: check your FXML file 'MeanOptions.fxml'.";

		/**
		 * MEANT TAB SET UP
		 */
		// Mean Tab Drag and Drop
		meanTabPane.setOnDragOver(new EventHandler<DragEvent>() {
			@Override
			public void handle(DragEvent event) {
				Dragboard db = event.getDragboard();
				if (db.hasFiles()) {
					event.acceptTransferModes(TransferMode.COPY);
				} else {
					event.consume();
				}
			}
		});
		// Dropping over surface
		meanTabPane.setOnDragDropped(new EventHandler<DragEvent>() {
			@Override
			public void handle(DragEvent event) {
				Dragboard db = event.getDragboard();
				boolean success = false;
				if (db.hasFiles()) {
					success = true;
					for (File file : db.getFiles()) {
						if (!meanfiles.contains(file))
							meanfiles.add(file);
					}
					mean_filenameObv.clear();
					for (File file : meanfiles)
						mean_filenameObv.add(file.getName());
				}
				event.setDropCompleted(success);
				event.consume();
			}
		});
		// Set List View with corresponding file names.
		mean_file_list.setItems(mean_filenameObv);
		// Allow User to remove file by clicking it.
		mean_file_list.getSelectionModel().selectedItemProperty()
				.addListener(new ChangeListener<String>() {
					@Override
					public void changed(ObservableValue<? extends String> arg0,
							String oldVal, String newVal) {
						mean_filenameObv.remove(newVal);
						mean_file_list.setItems(null);
						System.out.println(mean_filenameObv);
						System.out.println(meanfiles);
						mean_file_list.setItems(mean_filenameObv);
					}
				});

		// Set up Radio Toggles
		ToggleGroup radios = new ToggleGroup();
		rand_radio.setToggleGroup(radios);
		rand_radio.setSelected(true);
		randperm_radio.setToggleGroup(radios);
		// Set up Slider
		iterations_slide.setMin(0);
		iterations_slide.setMax(100000);
		iterations_max.setText("100000");
		// Update iterations label on sliding action.
		iterations_slide.valueProperty().addListener(
				new ChangeListener<Number>() {
					@Override
					public void changed(ObservableValue<? extends Number> ov,
							Number old_val, Number new_val) {
						iterations = new_val.intValue();
						iterations_min.setText(String.valueOf(new_val
								.intValue()));
					}
				});

		// Mean Submit Button Action
		sturmSubmit_button.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				// Begin Action
				meanSummary_label.setText("Working on Mean");
				String fullCommand = "";
				// Do for each File inputed.
				for (File inFile : meanfiles) {
					// Create Log for this file Iteration
					log = new File(inFile.getName().replace(".txt", "")
							+ "_STURM-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}
					// Set Up
					String eFactor = eFactor_field.getText();
					String eVal = eValue_field.getText();
					String interval = interval_field.getText();
					String cauchy = cauchy_field.getText();
					options = "";
					// Build Options from Radio Buttons
					if (rand_radio.isSelected())
						options += "-a/random/";
					else
						options += "-a/rand_perm/";
					// Build Epsilon
					if (eVal.length() < 1 && !(eFactor.length() < 1))
						options += "-f/" + eFactor + "/";
					if (!(eVal.length() < 1) && eFactor.length() < 1)
						options += "-e/" + eVal + "/";
					// Build Interval
					if (!(interval.length() < 1))
						options += "-i/" + interval + "/";
					// Always Build Iterations
					options += "-n/" + (int) iterations_slide.getValue() + "/";
					// Build Rooted, Permutate, Cauchy
					if (!meanRooted_toggle.isSelected())
						options += "-u/";
					if (perm_toggle.isSelected())
						options += "-p/";
					if (!(cauchy.length() < 1))
						options += "-c/" + cauchy + "/";
					// Retrieve output filename if exists otherwise create
					// default.
					String outFile = meanOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_MeanrawOut.txt/";

					// Build/Clean Call
					command = options + inFile.getPath();
					command = command.trim();
					// Call class' Main directly, log all System.out
					try {
						System.setOut(new PrintStream(log));
					} catch (FileNotFoundException e1) {
						System.out
								.println("Log message to show that the logger isn't logging the logs. "
										+ e1.toString());
					}
					try {
						CentroidMain.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute Mean "
								+ e1.toString());
					}
					// Update Progress Label
					fullCommand += command;
					meanSummary_label.setText(fullCommand);
				}
				// Open all generated files.
				fullCommand += "\nOpening Files...";
				meanSummary_label.setText(fullCommand);
				File folder = new File(System.getProperty("user.dir"));
				for (File f : folder.listFiles())
					if (f.getName().endsWith(".txt"))
						try {
							Desktop.getDesktop().open(f);
						} catch (IOException e1) {
							System.out.println("Couldn't open " + f
									+ e1.toString());
						}
				// Signal Mean has finished.
				fullCommand += "\nFiles Opened.";
				meanSummary_label.setText(fullCommand);
			}
		});
	}
}
