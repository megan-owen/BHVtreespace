package controllers;

import java.awt.Desktop;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.ResourceBundle;

import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.scene.Cursor;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.chart.ScatterChart;
import javafx.scene.control.Button;
import javafx.scene.control.ListView;
import javafx.scene.control.SplitPane;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleButton;
import javafx.scene.input.DragEvent;
import javafx.scene.input.Dragboard;
import javafx.scene.input.TransferMode;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import mdsj.IO;
import mdsj.MDSJ;

import org.python.util.PythonInterpreter;

import polyAlg.PolyMain;

public class GTPController {

	String options;
	String inFile;
	String outFile;
	String command;
	File folder;
	ArrayList<File> gtpfiles = new ArrayList<File>();
	ObservableList<String> gtp_filenameObv = FXCollections
			.observableArrayList("Input Files Here");
	File log;

	@FXML
	// ResourceBundle that was given to the FXMLLoader
	private ResourceBundle resources;

	@FXML
	// URL location of the FXML file that was given to the FXMLLoader
	private URL location;

	@FXML
	// fx:id="gtpTabPane"
	private SplitPane gtpTabPane; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtpOout_label"
	private Text gtpOout_label; // Value injected by FXMLLoader

	@FXML
	// fx:id="normalToggle"
	private ToggleButton normalToggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtpRooted_toggle"
	private ToggleButton gtpRooted_toggle; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtpOut_field"
	private TextField gtpOut_field; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtp_file_list"
	private ListView<String> gtp_file_list; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtpSubmit_button"
	private Button gtpSubmit_button; // Value injected by FXMLLoader

	@FXML
	// fx:id="gtpSummary_label"
	private Text gtpSummary_label; // Value injected by FXMLLoader

	@FXML
	// This method is called by the FXMLLoader when initialization is complete
	void initialize() {
		assert gtpTabPane != null : "fx:id=\"gtpTabPane\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtpOout_label != null : "fx:id=\"gtpOout_label\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert normalToggle != null : "fx:id=\"normalToggle\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtpRooted_toggle != null : "fx:id=\"gtpRooted_toggle\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtpOut_field != null : "fx:id=\"gtpOut_field\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtp_file_list != null : "fx:id=\"gtp_file_list\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtpSubmit_button != null : "fx:id=\"gtpSubmit_button\" was not injected: check your FXML file 'GTPOptions.fxml'.";
		assert gtpSummary_label != null : "fx:id=\"gtpSummary_label\" was not injected: check your FXML file 'GTPOptions.fxml'.";

		/**
		 * GTP TAB SET UP
		 */
		// GTP Tab Drag and Drop
		gtpTabPane.setOnDragOver(new EventHandler<DragEvent>() {
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
		gtpTabPane.setOnDragDropped(new EventHandler<DragEvent>() {
			@Override
			public void handle(DragEvent event) {
				Dragboard db = event.getDragboard();
				boolean success = false;
				if (db.hasFiles()) {
					success = true;
					//Add each file to file list, and show just names in listview
					for (File file : db.getFiles()) {
						if (!gtpfiles.contains(file))
							gtpfiles.add(file);
					}
					gtp_filenameObv.clear();
					for (File file : gtpfiles)
						gtp_filenameObv.add(file.getName());
				}
				event.setDropCompleted(success);
				event.consume();
			}
		});
		gtp_file_list.setItems(gtp_filenameObv);
		// Allow User to remove file by clicking it.
		gtp_file_list.getSelectionModel().selectedItemProperty()
		.addListener(new ChangeListener<String>() {
			@Override
			public void changed(ObservableValue<? extends String> arg0,
					String oldVal, String newVal) {
				gtp_filenameObv.remove(newVal);
				for (File file : gtpfiles)
					if (file.getName().equals(newVal))
						gtpfiles.remove(file);
				gtp_file_list.setItems(null);
				gtp_file_list.setItems(gtp_filenameObv);
			}
		});

		// GTP Submit Button Action
		gtpSubmit_button.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				String fullCommand = "";
				gtpSummary_label.setText("GTP Started...");
				for (File inFile : gtpfiles) {
					//Create new LOG file and set stream to print to it.
					log = new File(inFile.getName().replace(".txt", "")
							+ "_GTP-LOG.txt");
					PrintStream orig = System.out;
					try {
						PrintStream ps = new PrintStream(log);
						System.setOut(ps);
					} catch (FileNotFoundException e2) {
						System.out.println("Looks like that didn't work. "
								+ e2.toString());
					}

					options = "";
					// Build Options from Toggles
					if (normalToggle.isSelected())
						options += "-n/";
					if (!gtpRooted_toggle.isSelected())
						options += "-u/";
					//Additional toggles, removed.
					/**
					 * if(doubleToggle.isSelected()) options += "-d/";
					 * if(verboseToggle.isSelected()) options += "-v/";
					 **/
					// Retrieve output filename if exists.
					String outFile = gtpOut_field.getText();
					if (outFile.length() < 1)
						outFile = "default";
					if (outFile.endsWith(".txt"))
						outFile.replace(".txt", "");
					options += "-o/" + outFile + "_"
							+ inFile.getName().replace(".txt", "")
							+ "_GTPrawOut.txt/";
					command = options + inFile.getPath();
					command = command.trim();
					fullCommand += command + "\n";
					//Notify start of call
					System.out.println(fullCommand);
					gtpSummary_label.setText("Working on "
							+ fullCommand.replace("/", " "));
					//Send command directly to class Main
					try {
						PolyMain.main(command.split("/"));
					} catch (Exception e1) {
						System.out.println("Coulnd't execute GTP "
								+ e1.toString());
					}
					//Notify end of call
					gtpSummary_label.setText("Finished "
							+ fullCommand.replace("/", " "));
				}
				//Get current directory, used for managing generated files.
				File folder = new File(System.getProperty("user.dir"));
				//Run modded helper. Runs within Java to reduce dependencies.
				fullCommand += "\nCreating Matrix from GTP...";
				gtpSummary_label.setText(fullCommand);
				PythonInterpreter python = new PythonInterpreter();
				python.execfile("lib/modded_helper.py");
				//Runs MDS. Send all generated matrices to MDS lib.
				fullCommand += "\nRunning MDS from Matrix...";
				gtpSummary_label.setText(fullCommand);
				runMDS(folder);
				//Open all generated files, removed.
				/**
				fullCommand += "\nOpening all generated Files...";
				gtpSummary_label.setText(fullCommand);
					for (File f : folder.listFiles())
						if (f.getName().endsWith(".txt"))
							try {
								Desktop.getDesktop().open(f);
							} catch (IOException e1) {
								System.out.println("Couldn't open " + f
										+ e1.toString());
							}
				fullCommand += "\nFiles Opened.";
				 **/

				// Signal GTP Finish.
				gtpSummary_label.setText(fullCommand.replace("/", " "));
			}
		});
		/** GTP TAB SET UP END **/
	}
	//Helper method that receives the folder the program is currently running and searches through it for matrices.
	//Found matrices are sent to an MDS package which return the coordinates for each tree.
	//Temporarily these coordinates will be sent directly to a scatter plot in a new window.
	private static void runMDS(File folder){
		double[][] coords = null;
		double[][] matrix = null;
		for (File f : folder.listFiles())
			if (f.getName().endsWith("symMatrix.txt"))
				try {
					//Gather coords from lib.
					matrix = IO.read(f.getName());
					coords = MDSJ.classicalScaling(matrix);
					plotScatter(coords, f.getName().replace("_symMatrix.txt", ""));
				} catch (Exception e2) {
					e2.printStackTrace();
				}	
	}
	/**
	 * Receive coordinates and name of plot for that file.
	 * Creates a pop up window with a scatter chart of that files data.
	 * @param coords
	 * @param name
	 */
	private static void plotScatter(double[][] coords, String name){
		//Load scatter chart FXML and create a new Controller instance.
		ScatterChart scatter;
		ScatterController controller = new ScatterController(coords, name);
		try {
			scatter = FXMLLoader.load(ScatterController.class.getClass().getResource("/resources/scatterFrame.fxml"));
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		//Create new scene and set to a new pop up window.
		scatter = controller.getScatterChart();
		final Scene scene = new Scene(scatter);
		Stage popUp = new Stage();
		//Change cursor to crosshair- Visual cue to show points are to be clicked.
		scene.setOnMouseEntered(new EventHandler<Event>() {
			@Override
			public void handle(Event event) {
				scene.setCursor(Cursor.CROSSHAIR);
			}
		});
		scene.setOnMouseEntered(new EventHandler<Event>() {
			@Override
			public void handle(Event event) {
				scene.setCursor(Cursor.CROSSHAIR);
			}
		});
		popUp.setScene(scene);
		// CREATE AND SHOW WINDOW
		popUp.sizeToScene();
		popUp.setTitle(controller.chartName + " BHV Plot");
		popUp.show();
	}
}
