#ifndef ATHENA_READ_HPP_
#define ATHENA_READ_HPP_
#include <hdf5.h>
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <string>
#include "BootesArray.hpp"

using namespace std;
using namespace H5;
//https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html

/** Define variable index in the .athdf file
// -1 = don't read **/
const int ATHENA_VarIND_rho  = 0;
const int ATHENA_VarIND_pres = -1;
const int ATHENA_VarIND_vel1 = 1;
const int ATHENA_VarIND_vel2 = 2;
const int ATHENA_VarIND_vel3 = 3;
const int ATHENA_VarIND_Phi  = -1;
const int ATHENA_VarIND_Er   = -1;
const int ATHENA_VarIND_Fr1  = -1;
const int ATHENA_VarIND_Fr2  = -1;
const int ATHENA_VarIND_Fr3  = -1;
const int ATHENA_VarIND_Pr11 = -1;
const int ATHENA_VarIND_Pr22 = -1;
const int ATHENA_VarIND_Pr33 = -1;
const int ATHENA_VarIND_Pr12 = -1;
const int ATHENA_VarIND_Pr13 = -1;
const int ATHENA_VarIND_Pr23 = -1;
const int ATHENA_VarIND_Pr21 = -1;
const int ATHENA_VarIND_Pr31 = -1;
const int ATHENA_VarIND_Pr32 = -1;
const int ATHENA_VarIND_Er0  = -1;
const int ATHENA_VarIND_Fr01 = -1;
const int ATHENA_VarIND_Fr02 = -1;
const int ATHENA_VarIND_Fr03 = -1;
const int ATHENA_VarIND_sigma_s = -1;
const int ATHENA_VarIND_sigma_a = -1;
const int ATHENA_VarIND_sigma_p = -1;
const int ATHENA_VarIND_vp1 = -1;
const int ATHENA_VarIND_vp2 = -1;
const int ATHENA_VarIND_vp3 = -1;
const int ATHENA_VarIND_rhop = -1;


class OutputAthdf{
    public:
        OutputAthdf(string FILE_NAME){
            FILE_NAME_ = FILE_NAME;
            //setup();
            file = H5File( FILE_NAME_, H5F_ACC_RDONLY );
            ReadAttributes();          // Read necessary variables from Attributes
            ReadLogicLocations();      // Read the locations of each mesh block from DataSet
        }

        H5File file;
        unsigned int NumMeshBlocks;
        unsigned int NumCycle;
        unsigned int NumVariables;
        float Time;
        float RootGridX1[3];         // min, max and geometric ratio in X1
        float RootGridX2[3];         // min, max and geometric ratio in X2
        float RootGridX3[3];         // min, max and geometric ratio in X3
        unsigned int RootGridSize[3];
        unsigned int MeshBlockSize[3];
        BootesArray<int> LogicLocations;

        // Values of interface locations along x1/x2/x3-direction
        BootesArray<float> x1v;
        BootesArray<float> x2v;
        BootesArray<float> x3v;

        // Values of cell centers along x1/x2/x3-direction
        BootesArray<float> x1f;
        BootesArray<float> x2f;
        BootesArray<float> x3f;

        BootesArray<float> rho;
        BootesArray<float> temp;
        BootesArray<float> pres;
        BootesArray<float> Phi;
        BootesArray<float> vel1;
        BootesArray<float> vel2;
        BootesArray<float> vel3;

    void ByteSwap32Int(unsigned int *value, int length);
    void ByteSwap32Int(unsigned int &value);
    void ReadAttributes();
    void SetupCartesianGrid(float lscale);
    void ReadLogicLocations();
    void ReadEachMeshBlock(string DataSetName, unsigned int hdf5Start[5], unsigned int hdf5Select[5], unsigned int outputshape[3], BootesArray<float> &data_out);
    void ReadAndCombineMeshBlocks(unsigned int ATHENA_VarIND, BootesArray<float> &data, float ScaleVal);
    void DoRead(int VarIND, BootesArray<float> &arr, string VarName, float ScaleVal);
    void FillValue(BootesArray<float> &arr, float fill_value);
    template <typename T>
    void getAttribute(string att_name, T &value);
    void setup();
    void CalcTemperature();
    void updatefile(string new_fn);
    void close();

    private:
        string FILE_NAME_;
};


void getAttributeString(H5File file, string att_name, H5std_string &value_str, const int NumVars){
    Attribute att = file.openAttribute(att_name);
    DataType type = att.getDataType();

    att.read(type, &value_str);
}


void OutputAthdf::ByteSwap32Int(unsigned int &value){
    value = __builtin_bswap32 (value);
}


void OutputAthdf::ByteSwap32Int(unsigned int *value, int length){
    for (int ii = 0; ii < length; ii ++ ){
        value[ii] = __builtin_bswap32 (value[ii]);
    }
}


template <typename T>
void OutputAthdf::getAttribute(string att_name, T &value){
    Attribute att = file.openAttribute(att_name);
    DataType type = att.getDataType();
    att.read(type, &value);
}


void OutputAthdf::ReadAttributes(){
        getAttribute<unsigned int>("NumCycles", NumCycle);
        ByteSwap32Int(NumCycle);
        getAttribute<unsigned int>("NumMeshBlocks", NumMeshBlocks);
        ByteSwap32Int(NumMeshBlocks);
        getAttribute<unsigned int>("NumVariables", NumVariables);
        ByteSwap32Int(NumVariables);
        double Time_temp;
        double RootGridX1_temp[3];
        double RootGridX2_temp[3];
        double RootGridX3_temp[3];
        getAttribute<double>("Time", Time_temp);
        getAttribute<double>("RootGridX1", *RootGridX1_temp);
        getAttribute<double>("RootGridX2", *RootGridX2_temp);
        getAttribute<double>("RootGridX3", *RootGridX3_temp);
        Time = (float) Time_temp;
        for (int ii = 0; ii < 3; ii ++ ){
            RootGridX1[ii] = (float) RootGridX1_temp[ii];
            RootGridX2[ii] = (float) RootGridX2_temp[ii];
            RootGridX3[ii] = (float) RootGridX3_temp[ii];
        }
        getAttribute<unsigned int>("RootGridSize", *RootGridSize);
        ByteSwap32Int(RootGridSize, 3);
        getAttribute<unsigned int>("MeshBlockSize", *MeshBlockSize);
        ByteSwap32Int(MeshBlockSize, 3);
    }


void OutputAthdf::SetupCartesianGrid(float lscale){
    x1f.NewBootesArray(RootGridSize[0] + 1);
    x2f.NewBootesArray(RootGridSize[1] + 1);
    x3f.NewBootesArray(RootGridSize[2] + 1);
    x1v.NewBootesArray(RootGridSize[0]);
    x2v.NewBootesArray(RootGridSize[1]);
    x3v.NewBootesArray(RootGridSize[2]);

    RootGridX1[1] *= lscale;
    RootGridX2[1] *= lscale;
    RootGridX3[1] *= lscale;
    RootGridX1[0] *= lscale;
    RootGridX2[0] *= lscale;
    RootGridX3[0] *= lscale;
    float dx1 = (RootGridX1[1] - RootGridX1[0]) / RootGridSize[0];
    float dx2 = (RootGridX2[1] - RootGridX2[0]) / RootGridSize[1];
    float dx3 = (RootGridX3[1] - RootGridX3[0]) / RootGridSize[2];

    for (int ii = 0; ii < RootGridSize[0] + 1; ii ++){
        x1f(ii) = (RootGridX1[0] + dx1 * ii);
    }
    for (int ii = 0; ii < RootGridSize[1] + 1; ii ++){
        x2f(ii) = (RootGridX2[0] + dx2 * ii);
    }
    for (int ii = 0; ii < RootGridSize[2] + 1; ii ++){
        x3f(ii) = (RootGridX3[0] + dx3 * ii);
    }

    for (int ii = 0; ii < RootGridSize[0]; ii ++){
        x1v(ii) = (x1f(ii) + x1f(ii + 1)) / 2.;
    }
    for (int ii = 0; ii < RootGridSize[1]; ii ++){
        x2v(ii) = (x2f(ii) + x2f(ii + 1)) / 2.;
    }
    for (int ii = 0; ii < RootGridSize[2]; ii ++){
        x3v(ii) = (x3f(ii) + x3f(ii + 1)) / 2.;
    }
}


void OutputAthdf::ReadLogicLocations(){
    // >>> GET the locations of each mesh block >>>
    DataSet dataset_logical_loc = file.openDataSet("LogicalLocations");
    H5T_class_t type_class_logical_loc = dataset_logical_loc.getTypeClass();
    IntType intype = dataset_logical_loc.getIntType();
    /*
    * Get dataspace of the dataset.
    */
    DataSpace dataspace_logical_loc = dataset_logical_loc.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    int rank_logical_loc = dataspace_logical_loc.getSimpleExtentNdims();
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t d_logical_loc;
    hsize_t dims_out_logical_loc[2];
    int ndims_logical_loc = dataspace_logical_loc.getSimpleExtentDims( dims_out_logical_loc, NULL);

    /*
    * Define hyperslab in the dataset; implicitly giving strike and
    * block NULL.
    */
    hsize_t      offset_logical_loc[2];   // hyperslab offset in the file
    hsize_t      count_logical_loc[2];    // size of the hyperslab in the file
    offset_logical_loc[0] = 0;               // starting location of hyperslab in 1st axis
    offset_logical_loc[1] = 0;               // starting location of hyperslab in 2nd axis
    count_logical_loc[0]  = NumMeshBlocks;   // number of locations to be selected in 1st axis
    count_logical_loc[1]  = 3;               // number of locations to be selected in 2nd axis
    dataspace_logical_loc.selectHyperslab( H5S_SELECT_SET, count_logical_loc, offset_logical_loc );
    /*
    * Define the memory dataspace, used for output
    */
    hsize_t     dimsm_logical_loc[2];              /* memory space dimensions */
    dimsm_logical_loc[0] = NumMeshBlocks;
    dimsm_logical_loc[1] = 3;
    DataSpace memspace_logical_loc( 2, dimsm_logical_loc );
    /*
    * Define memory hyperslabm used for output
    */
    hsize_t      offset_out_logical_loc[2];   // hyperslab offset in memory
    hsize_t      count_out_logical_loc[2];    // size of the hyperslab in memory
    offset_out_logical_loc[0] = 0;
    offset_out_logical_loc[1] = 0;
    count_out_logical_loc[0]  = NumMeshBlocks;
    count_out_logical_loc[1]  = 3;
    memspace_logical_loc.selectHyperslab( H5S_SELECT_SET, count_out_logical_loc, offset_out_logical_loc );

    /*
    * Read data from hyperslab in the file into the hyperslab in
    * memory and display the data.
    */
    int *LogicLocationsTemp = new int[NumMeshBlocks * 3];  // shape1 = 3
    dataset_logical_loc.read( LogicLocationsTemp, PredType::NATIVE_INT, memspace_logical_loc, dataspace_logical_loc );

    LogicLocations.NewBootesArray(NumMeshBlocks, 3);
    for (int jj = 0; jj < NumMeshBlocks; jj++){
        for (int ii = 0; ii < 3; ii++){
            LogicLocations(jj, ii) = LogicLocationsTemp[jj * 3 + ii];        // correct iteration;
        }
    }
    delete []LogicLocationsTemp;
}


void OutputAthdf::ReadEachMeshBlock(string DataSetName, unsigned int hdf5Start[5], unsigned int hdf5Select[5], unsigned int outputshape[3], BootesArray<float> &data_out){
    /**
        Read the selected 3D data cube from the 5-dimensional datacube "Prim" in .athdf file

        indend in:  H5File file
        indend in:  string DataSetName
        indend in:  unsigned int hdf5Start[5]
        indend in:  unsigned int hdf5Select[5]
        indend out: BootesArray<float> data_out
    **/
    DataSet dataset = file.openDataSet(DataSetName);
    H5T_class_t type_class = dataset.getTypeClass();

    FloatType floatype = dataset.getFloatType();
    /*
    * Get order of datatype and print message if it's a little endian.
    */
    H5std_string order_string;
    H5T_order_t order = floatype.getOrder( order_string );

    size_t size = floatype.getSize();

    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */

    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[5];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

    /*
    * Define hyperslab in the dataset; implicitly giving strike and
    * block NULL.
    */
    hsize_t      offset[5];      // hyperslab offset in the file
    hsize_t      count[5];       // size of the hyperslab in the file
    offset[0] = hdf5Start[0];               // starting location of hyperslab in 1st axis
    offset[1] = hdf5Start[1];               // starting location of hyperslab in 2nd axis
    offset[2] = hdf5Start[2];               // starting location of hyperslab in 3rd axis
    offset[3] = hdf5Start[3];               // starting location of hyperslab in 4th axis
    offset[4] = hdf5Start[4];               // starting location of hyperslab in 5th axis
    count[0]  = hdf5Select[0];               // number of locations to be selected in 1st axis
    count[1]  = hdf5Select[1];               // number of locations to be selected in 2nd axis
    count[2]  = hdf5Select[2];             // number of locations to be selected in 3rd axis
    count[3]  = hdf5Select[3];             // number of locations to be selected in 4th axis
    count[4]  = hdf5Select[4];             // number of locations to be selected in 5th axis
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset);

    /*
    * Define the memory dataspace.
    */
    hsize_t dimsm[3];
    dimsm[0] = outputshape[0];
    dimsm[1] = outputshape[1];
    dimsm[2] = outputshape[2];
    unsigned int RANK_OUT = 3;
    DataSpace memspace( RANK_OUT, dimsm );

    /*
    * Define memory hyperslab.
    */
    hsize_t      offset_out[3];   // hyperslab offset in memory
    hsize_t      count_out[3];    // size of the hyperslab in memory
    offset_out[0] = 0;
    offset_out[1] = 0;
    offset_out[2] = 0;
    count_out[0]  = outputshape[0];
    count_out[1]  = outputshape[1];
    count_out[2]  = outputshape[2];
    memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );
    float *data_out_temp = new float[outputshape[0] * outputshape[1] * outputshape[2]];
    dataset.read( data_out_temp, PredType::NATIVE_FLOAT, memspace, dataspace );
    //#pragma omp parallel for schedule(static)
    for (int kk = 0; kk < outputshape[2]; kk ++){
        for (int jj = 0; jj < outputshape[1]; jj ++){
            for (int ii = 0; ii < outputshape[0]; ii ++){
               data_out(kk, jj, ii) = data_out_temp[(kk * outputshape[1] + jj) * outputshape[0] + ii];
               // data_out(kk, jj, ii) = data_out_temp[kk][jj][ii];   x1 + size1_*(x2 + size2_*x3)
            }
        }
    }
    //#pragma omp barrier
    delete []data_out_temp;

}


void OutputAthdf::ReadAndCombineMeshBlocks(unsigned int VarIND, BootesArray<float> &data, float ScaleVal){
    for (unsigned int MB_IND = 0; MB_IND < NumMeshBlocks; MB_IND ++){
        // read out data from file
        unsigned int primstart[5]  = {VarIND, MB_IND, 0, 0, 0};
        unsigned int primselect[5] = {1, 1, MeshBlockSize[2], MeshBlockSize[1], MeshBlockSize[0]};
        unsigned int outputshape[3] = {MeshBlockSize[2], MeshBlockSize[1], MeshBlockSize[0]};

        BootesArray<float> data_out;
        data_out.NewBootesArray(MeshBlockSize[2], MeshBlockSize[1], MeshBlockSize[0]);
        ReadEachMeshBlock("prim", primstart, primselect, outputshape, data_out);
        //assign the data to their location(s) on mesh
        int MB_loc2 = LogicLocations(MB_IND, 2) * MeshBlockSize[2];
        int MB_loc1 = LogicLocations(MB_IND, 1) * MeshBlockSize[1];
        int MB_loc0 = LogicLocations(MB_IND, 0) * MeshBlockSize[0];

        // #pragma omp parallel for schedule(static)
        for (int kk = 0; kk < MeshBlockSize[2]; kk ++){
            for (int jj = 0; jj < MeshBlockSize[1]; jj ++){
                for (int ii = 0; ii < MeshBlockSize[0]; ii ++){
                    data(MB_loc2 + kk, MB_loc1 + jj, MB_loc0 + ii) = data_out(kk, jj, ii) * ScaleVal;
                }
            }
        }
        //#pragma omp barrier
    }
}


void OutputAthdf::DoRead(int VarIND, BootesArray<float> &arr, string VarName, float ScaleVal){
    /* read out density */        // Index order (r, theta, phi)
    if (VarIND  != -1) {
        arr.NewBootesArray(RootGridSize[2], RootGridSize[1], RootGridSize[0]);
        ReadAndCombineMeshBlocks(VarIND, arr, ScaleVal);
    }
    else {
        cout << "No " << VarName << " Data" << endl << flush;
    }
}


void OutputAthdf::FillValue(BootesArray<float> &arr, float fill_value){
    arr.NewBootesArray(RootGridSize[2], RootGridSize[1], RootGridSize[0]);
    for (int kk = 0; kk < RootGridSize[2]; kk ++){
        for (int jj = 0; jj < RootGridSize[1]; jj ++){
            for (int ii = 0; ii < RootGridSize[0]; ii ++){
                arr(kk, jj, ii) = fill_value;
            }
        }
    }
}


void OutputAthdf::CalcTemperature(){
    temp.NewBootesArray(RootGridSize[2], RootGridSize[1], RootGridSize[0]);
    for (int kk = 0; kk < RootGridSize[2]; kk ++){
        for (int jj = 0; jj < RootGridSize[1]; jj ++){
            for (int ii = 0; ii < RootGridSize[0]; ii ++){
                temp(kk, jj, ii) = pres(kk, jj, ii) * hydro_mu * mH / (rho(kk, jj, ii) * kb);
            }
        }
    }
}


void OutputAthdf::updatefile(string new_fn){
    //cout << "warning: " << FILE_NAME_ << " is being replaced by " << new_fn << endl << flush;
    file.close();
    /** re-initialize **/
    FILE_NAME_ = new_fn;
    //setup();
    file = H5File( FILE_NAME_, H5F_ACC_RDONLY );
    // clean up old file
    LogicLocations.clean();
    x1v.clean();
    x2v.clean();
    x3v.clean();
    x1f.clean();
    x2f.clean();
    x3f.clean();
    rho.clean();
    pres.clean();
    temp.clean();
    Phi.clean();
    vel1.clean();
    vel2.clean();
    vel3.clean();
    // read in new ones
    ReadAttributes();          // Read necessary variables from Attributes
    ReadLogicLocations();      // Read the locations of each mesh block from DataSet
}


void OutputAthdf::close(){
    file.close();
}
#endif // ATHENA_READ_HPP_
