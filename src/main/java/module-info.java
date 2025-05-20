module ca.smithlab.vandalovejoy.ecsfinder {
    requires javafx.controls;
    requires javafx.fxml;
            
                            
    opens ca.smithlab.vandalovejoy.ecsfinder to javafx.fxml;
    exports ca.smithlab.vandalovejoy.ecsfinder;
}