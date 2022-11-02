//
// Created by joeho on 13/08/2022.
//

#include "ioparser.h"
using namespace std;

IOParser::IOParser(string out_disp_fpath, string out_log_fpath) {
    /**
     * @brief instantiate the IOParser object.
     * @param setting_fpath the json file to open and read settings from.
     * @param out_fpath the binary file fp to write simulation output to.
     * 
     */
    io_out_disp_fpath = out_disp_fpath;
    io_out_log_fpath = out_log_fpath;
}

void::IOParser::OpenOutFile() {
    /**
     * @brief Opens the output binary file stream.
     */
    io_out_disp_file.open(io_out_disp_fpath, ios::out | ios::binary);
}

void::IOParser::WriteStateToOutFile(std::list<particle_type*> particle_list, int n) {
    /**
     * @brief Writes the given state to the IOParser binary out file, as an
     * array of x coords, y coords and then colours, all as floats. Does
     * not open or close the filestream, simply flushes the output stream.
     */

    // create temporary arrays for writing
    float* arr_x = new float[n];
    float* arr_y = new float[n];
    float* arr_colour = new float[n];

    // fill temp arrays from particle structs in state
    int i = 0;
    for (particle_type* &p : particle_list) {
        arr_x[i] = p->x;
        arr_y[i] = p->y;
        arr_colour[i] = (float)p->colour;
        i++;
    }

    float nf = (float)n;

    // write number or particles for this state, write coords, and colours
    io_out_disp_file.write((char *)&nf, sizeof(float));
    io_out_disp_file.write((char *)arr_x, sizeof(float) * n);
    io_out_disp_file.write((char *)arr_y, sizeof(float) * n);
    io_out_disp_file.write((char *)arr_colour, sizeof(float) * n);

    // flush output stream, but keep file open
    io_out_disp_file.flush();

    // clean temp arrays
    delete[] arr_x;
    delete[] arr_y;
    delete[] arr_colour;
}

void::IOParser::CloseOutFile() {
    /**
     * @brief Closes the output binary file stream.
     */
    io_out_disp_file.close();
}




