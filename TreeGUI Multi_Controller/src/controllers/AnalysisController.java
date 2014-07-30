package controllers;

import java.awt.Desktop;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.ResourceBundle;

import distances.Analysis;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.scene.control.ListView;
import javafx.scene.control.Slider;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleButton;
import javafx.scene.control.ToggleGroup;
import javafx.scene.input.DragEvent;
import javafx.scene.input.Dragboard;
import javafx.scene.input.TransferMode;
import javafx.scene.text.Text;

public class AnalysisController {
	String options;
	String inFile;
	String outFile;
	String command;
	File folder;
	ArrayList<File> analysisfiles = new ArrayList<File>();
	ObservableList<String> analysis_filenameObv = FXCollections
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
	// fx:id="sampleVal"
	private Text sampleVal; // Value injected by FXMLLoader

	@FXML
	// fx:id="analysisSummary_label"
	private Text analysisSummary_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="sample_slider"
	private Slider sample_slider; // Value injected by FXMLLoader

	@FXML
	// fx:id="analysisOut_field"
	private TextField analysisOut_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="sampleMax"
	private Text sampleMax; // Value injected by FXMLLoader

	@FXML
	// fx:id="analysisTabPane"
	private SplitPane analysisTabPane; // Value injected by FXMLLoader

	@FXML
	// fx:id="length_toggle"
	private ToggleButton length_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="analysis_file_list"
	private ListView<String> analysis_file_list; // Value injected by FXMLLoader

	@FXML
	// fx:id="sample_toggle"
	private ToggleButton sample_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="project_toggle"
	private ToggleButton project_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="analysisOut_label"
	private Text analysisOut_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="sampleMin"
	private Text sampleMin; // Value injected by FXMLLoader

	@FXML
	// fx:id="geo_toggle"
	private ToggleButton geo_toggle; // Value injected by FXMLLoader

	@FXML
	// This method is called by the FXMLLoader when initialization is complete
	void initialize() {
		assert sampleVal != null : "fx:id=\"sampleVal\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert analysisSummary_label != null : "fx:id=\"analysisSummary_label\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert sample_slider != null : "fx:id=\"sample_slider\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert analysisOut_field != null : "fx:id=\"analysisOut_field\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert sampleMax != null : "fx:id=\"sampleMax\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert analysisTabPane != null : "fx:id=\"analysisTabPane\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert length_toggle != null : "fx:id=\"length_toggle\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert analysis_file_list != null : "fx:id=\"analysis_file_list\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert sample_toggle != null : "fx:id=\"sample_toggle\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert project_toggle != null : "fx:id=\"project_toggle\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert analysisOut_label != null : "fx:id=\"analysisOut_label\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert sampleMin != null : "fx:id=\"sampleMin\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";
		assert geo_toggle != null : "fx:id=\"geo_toggle\" was not injected: check your FXML file 'AnalysisOptions.fxml'.";

		/**
		 * ANALYSIS TAB BUILD
		 */
		// Drag and Drop
		analysisTabPane.setOnDragOver(new EventHandler<DragEvent>() {
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
		analysisTabPane.setOnDragDropped(new EventHandler<DragEvent>() {
			@Override
			public void handle(DragEvent event) {
				Dragboard db = event.getDragboard();
				boolean success = false;
				if (db.hasFiles()) {
					success = true;
					for (File file : db.getFiles()) {
						if (!analysisfiles.contains(file))
							analysisfiles.add(file);
					}
					analysis_filenameObv.clear();
					for (File file : analysisfiles)
						analysis_filenameObv.add(file.getName());
				}
				event.setDropCompleted(success);
				event.consume();
			}
		});
		analysis_file_list.setItems(analysis_filenameObv);
		// Allow User to remove file by clicking it.
		analysis_file_list.getSelectionModel().selectedItemProperty()
				.addListener(new ChangeListener<String>() {
					@Override
					public void changed(ObservableValue<? extends String> arg0,
							String oldVal, String newVal) {
						analysis_filenameObv.remove(newVal);
						analysis_file_list.setItems(null);
						System.out.println(analysis_filenameObv);
						System.out.println(analysisfiles);
						analysis_file_list.setItems(analysis_filenameObv);
					}
				});
		// Toggle Set up in Group
		ToggleGroup toggles = new ToggleGroup();
		geo_toggle.setToggleGroup(toggles);
		project_toggle.setToggleGroup(toggles);
		length_toggle.setToggleGroup(toggles);
		sample_toggle.setToggleGroup(toggles);
		// Set up slider action.
		sample_slider.valueProperty().addListener(new ChangeListener<Number>() {
			@Override
			public void changed(ObservableValue<? extends Number> ov,
					Number old_val, Number new_val) {
				sampleVal.setText(String.valueOf(new_val.doubleValue())
						.substring(0, 4));
			}
		});
		// Analysis Toggle Actions
		geo_toggle.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				String fullCommand = "";
				analysisSummary_label.setText("Dual GTP Started...");
				for (int c = 0; c < analysisfiles.size(); c += 2) {
					File inFile = analysisfiles.get(c);
					File secondFile = analysisfiles.get(c + 1);
					log = new File(inFile.getName().replace(".txt", "")
							+ "_ANALYSIS-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}
					options = "-a/gtp_twofiles";
					options += "/-f/" + secondFile.getPath();

					// Retrieve output filename if exists.
					String outFile = analysisOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "/-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_rawOut.txt/";
					// Build/Clean Call
					command = options + inFile.getPath();
					command = command.trim();
					fullCommand += command + "\n";
					// Call class Main directly and log all System.out
					try {
						System.setOut(new PrintStream(log));
					} catch (FileNotFoundException e1) {
						System.out
								.println("Log message to show that the logger isn't logging the logs."
										+ e1.toString());
					}
					try {
						Analysis.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute dual GTP "
								+ e1.toString());
					}
					analysisSummary_label.setText("Finished "
							+ fullCommand.replace("/", " "));
				}
				// Open all generated files.
				fullCommand += "\nOpening Files...";
				analysisSummary_label.setText(fullCommand);
				File folder = new File(System.getProperty("user.dir"));
				for (File f : folder.listFiles())
					if (f.getName().endsWith(".txt"))
						try {
							Desktop.getDesktop().open(f);
						} catch (IOException e1) {
							System.out.println("Couldn't open " + f
									+ e1.toString());
						}
				fullCommand += "\nFiles Opened.";
				analysisSummary_label.setText(fullCommand);
			}
		});
		project_toggle.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				String fullCommand = "";
				analysisSummary_label.setText("Project Started...");
				for (int c = 0; c < analysisfiles.size(); c += 2) {
					File inFile = analysisfiles.get(c);
					File secondFile = analysisfiles.get(c + 1);
					log = new File(inFile.getName().replace(".txt", "")
							+ "_ANALYSIS-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}
					options = "-a/project";
					options += "/-f/" + secondFile.getPath();

					// Retrieve output filename if exists.
					String outFile = analysisOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "/-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_AnalysisRawOut.txt/";
					// Build/Clean Call
					command = options + inFile.getPath();
					command = command.trim();
					fullCommand += command + "\n";
					// Call class Main directly and log all System.out
					try {
						System.setOut(new PrintStream(log));
					} catch (FileNotFoundException e1) {
						System.out
								.println("Log message to show that the logger isn't logging the logs."
										+ e1.toString());
					}
					try {
						Analysis.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute projection "
								+ e1.toString());
					}
					analysisSummary_label.setText("Finished "
							+ fullCommand.replace("/", " "));
				}
				// Open all generated files.
				fullCommand += "\nOpening Files...";
				analysisSummary_label.setText(fullCommand);
				File folder = new File(System.getProperty("user.dir"));
				for (File f : folder.listFiles())
					if (f.getName().endsWith(".txt"))
						try {
							Desktop.getDesktop().open(f);
						} catch (IOException e1) {
							System.out.println("Couldn't open " + f
									+ e1.toString());
						}
				fullCommand += "\nFiles Opened.";
				analysisSummary_label.setText(fullCommand);
			}
		});
		length_toggle.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				String fullCommand = "";
				analysisSummary_label.setText("Metrics Started...");
				for (File inFile : analysisfiles) {
					log = new File(inFile.getName().replace(".txt", "")
							+ "_ANALYSIS-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}
					options = "-a/lengths/-a/topology_count/-a/split_count";

					// Retrieve output filename if exists.
					String outFile = analysisOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "/-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_AnalysisRawOut.txt/";
					// Build/Clean Call
					command = options + inFile.getPath();
					command = command.trim();
					fullCommand += command + "\n";
					// Call class Main directly and log all System.out
					try {
						System.setOut(new PrintStream(log));
					} catch (FileNotFoundException e1) {
						System.out
								.println("Log message to show that the logger isn't logging the logs."
										+ e1.toString());
					}
					try {
						Analysis.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute metrics "
								+ e1.toString());
					}
					analysisSummary_label.setText("Finished "
							+ fullCommand.replace("/", " "));
				}
				// Open all generated files.
				fullCommand += "\nOpening Files...";
				analysisSummary_label.setText(fullCommand);
				File folder = new File(System.getProperty("user.dir"));
				for (File f : folder.listFiles())
					if (f.getName().endsWith(".txt"))
						try {
							Desktop.getDesktop().open(f);
						} catch (IOException e1) {
							System.out.println("Couldn't open " + f
									+ e1.toString());
						}
				fullCommand += "\nFiles Opened.";
				analysisSummary_label.setText(fullCommand);
			}
		});
		sample_toggle.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				String fullCommand = "";
				analysisSummary_label.setText("Sample Point Started...");
				for (File inFile : analysisfiles) {
					log = new File(inFile.getName().replace(".txt", "")
							+ "_ANALYSIS-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}
					options = "-a/sample_point/-s/" + sample_slider.getValue();

					// Retrieve output filename if exists.
					String outFile = analysisOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "/-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_AnalysisRawOut.txt/";
					// Build/Clean Call
					command = options + inFile.getPath();
					command = command.trim();
					fullCommand += command + "\n";
					// Call class Main directly and log all System.out
					try {
						System.setOut(new PrintStream(log));
					} catch (FileNotFoundException e1) {
						System.out
								.println("Log message to show that the logger isn't logging the logs."
										+ e1.toString());
					}
					try {
						Analysis.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute metrics "
								+ e1.toString());
					}
					analysisSummary_label.setText("Finished "
							+ fullCommand.replace("/", " "));
				}
				// Open all generated files.
				fullCommand += "\nOpening Files...";
				analysisSummary_label.setText(fullCommand);
				File folder = new File(System.getProperty("user.dir"));
				for (File f : folder.listFiles())
					if (f.getName().endsWith(".txt"))
						try {
							Desktop.getDesktop().open(f);
						} catch (IOException e1) {
							System.out.println("Couldn't open " + f
									+ e1.toString());
						}
				fullCommand += "\nFiles Opened.";
				analysisSummary_label.setText(fullCommand);
			}
		});
	}
}
