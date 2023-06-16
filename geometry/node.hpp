class Node {
    private:
        short ID;
        float x_coordinate;
        float y_coordinate;
        float z_coordinate;

    public:
        Node(short identifier, float x_value, float y_value, float z_value){
            ID = identifier;
            x_coordinate = x_value;
            y_coordinate = y_value;
            z_coordinate = z_value;
        }

        void set_ID(short identifier){
            ID = identifier;
        }
        short get_ID(){
            return ID;
        }

        void set_x_coordinate(float x_value){
            x_coordinate = x_value;
        }
        float get_x_coordinate(){
            return x_coordinate;
        }

        void set_y_coordinate(float y_value){
            y_coordinate = y_value;
        }
        float get_y_coordinate(){
            return y_coordinate;
        }

        void set_z_coordinate(float z_value){
            z_coordinate = z_value;
        }
        float get_z_coordinate(){
            return z_coordinate;
        }
};
