package controllers;
import javafx.event.Event;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.scene.Node;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.input.MouseEvent;

public class ScatterController {
	//Required FXML to denote this is a controller of a ScatterChart
	@FXML
	private ScatterChart<Number, Number> scatterChart;
	//Chart var
	public String chartName;
	private double[][] coords;
	public XYChart.Series points = new XYChart.Series();
	//Constructor receives a name and raw coords for the scatter chart.
	public ScatterController(double[][] coord, String name){
		chartName = name;
		this.coords = coord;
		plotScatter();
	}
	//Takes raw coords and creates a new series which is added as data to the
	//Scatterchart.
	private void plotScatter(){
		scatterChart = new ScatterChart<Number, Number>(getXaxis(),getYaxis());
		for (int i = 0; i < coords[0].length; i++) {
			points.getData().add(new XYChart.Data(coords[0][i],coords[1][i]));
		}
		applyOnClickData(points);
		scatterChart.getData().setAll(points);
	}
	private void applyOnClickData(final Series series) {
		Node node = series.getNode();
		node.setOnMouseClicked(new EventHandler<MouseEvent>(){
			@Override
			public void handle(MouseEvent event) {
				System.out.println(series.getData());
			}
			
		});
	}
	//Simple helper methods that take the raw coordinates and return Axis
	//objects that capture the complete data in the chart.
	private NumberAxis getXaxis(){
		NumberAxis xAxis = null;
		double smallest = coords[0][0];
		double largest = coords[0][0];
		for(int c = 0; c < coords[0].length; c++){
			if(coords[0][c] < smallest)
				smallest = coords[0][c];
			if(coords[0][c] > largest)
				largest = coords[0][c];
		}
		xAxis = new NumberAxis(smallest, largest, (largest-smallest)/coords[0].length);
		return xAxis;
	}
	private NumberAxis getYaxis(){
		NumberAxis yAxis = null;
		double smallest = coords[1][0];
		double largest = coords[1][0];
		for(int c = 0; c < coords[0].length; c++){
			if(coords[1][c] < smallest)
				smallest = coords[1][c];
			if(coords[1][c] > largest)
				largest = coords[1][c];
		}
		yAxis = new NumberAxis(smallest, largest, (largest-smallest)/coords[0].length);
		return yAxis;
	}
	//Simple Get Method
	public ScatterChart getScatterChart(){
		return scatterChart;
	}
}
