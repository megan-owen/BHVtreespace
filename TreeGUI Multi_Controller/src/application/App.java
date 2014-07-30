/**
 * Written by Marlon Figueroa
 * June2-, 2014
 */

package application;

import java.io.IOException;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.image.Image;
import javafx.stage.Stage;

public class App extends Application {
	// Called before initialize()
	// Loads FXML
	// Creates Window from FXML
	@Override
	public void start(Stage primaryWindow) {
		Parent root;
		try {
			root = FXMLLoader.load(getClass().getResource("/resources/guiFrame.fxml"));
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		Scene scene = new Scene(root);

		// CREATE AND SHOW WINDOW
		primaryWindow.setScene(scene);
		primaryWindow.sizeToScene();
		primaryWindow.setTitle("Owen Treespace");
		primaryWindow.setWidth(800);
		primaryWindow.setHeight(800);
		primaryWindow.setResizable(false);
		primaryWindow.show();
		primaryWindow.getIcons().add(new Image("/resources/icon.png"));
		// Prevent Unwanted Close event.
		// NOT NEEDED
		/**
		 * Platform.setImplicitExit(false); primaryWindow.setOnCloseRequest(new
		 * EventHandler<WindowEvent>() {
		 * 
		 * @Override public void handle(WindowEvent event) { event.consume(); }
		 *           });
		 **/
	}

	// Launch the App
	public static void main(String[] args) {
		launch(args);
	}
}
