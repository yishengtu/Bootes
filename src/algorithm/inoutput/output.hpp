#ifndef OUTPUT_HPP_
#define OUTPUT_HPP_

#include <vector>
#include <hdf5.h>
#include <H5Cpp.h>
#include <iostream>
#include <string>
#include "../BootesArray.hpp"

using namespace std;
using namespace H5;

class Output{
    public:
        Output(string fn, char mode){
            fn_ = fn;
            mode_ = mode;
            if (mode == 'w'){
                file_ = new H5File(fn, H5F_ACC_TRUNC);
            }
            else if (mode == 'r'){
                cout << "TODO-impliment reading output" << endl;
            }
            else{
                cout << "mode not reconized" << endl;
                throw 1;
            }
        }

    template<typename T>
    void writeattribute(T *attr, string name, PredType hdf5_type, unsigned int FSPACE_ATT){
        /**
            T* att: attribute to be written
            name: name of attribute
            hdf5_type: type in hdf5 file
            FSPACE ATT: assume attributes are always 1-dimensional, the length of the 1D array (length of att)
        **/

        unsigned int FSPACE_ATT_RANK = 1;
        hsize_t fdim_att[] = {FSPACE_ATT}; // dim sizes of ds (on disk)
        DataSpace fspace_att( FSPACE_ATT_RANK, fdim_att );

        Attribute* att = new Attribute(file_->createAttribute(name, hdf5_type, fspace_att));
        att->write( hdf5_type, attr);
        delete att;
    }

    template<typename T>
    void write1Ddataset(BootesArray<T> &datain, string name, PredType hdf5_type){
        if (datain.dimension() != 1){
            throw 1;
        }

        const unsigned int FSPACE_DIM1 = datain.shape()[0];
        unsigned int FSPACE_RANK = 1;

        hsize_t fdim[] = {FSPACE_DIM1}; // dim sizes of ds (on disk)
        DataSpace fspace( FSPACE_RANK, fdim );

        DataSet* dataset = new DataSet(file_->createDataSet(name, hdf5_type, fspace));
        /*
         * Select hyperslab for the dataset in the file, using 3x2 blocks,
         * (4,3) stride and (2,4) count starting at the position (0,1).
         */
        hsize_t start[1];  // Start of hyperslab
        hsize_t stride[1]; // Stride of hyperslab
        hsize_t count[1];  // Block count
        hsize_t block[1];  // Block sizes
        start[0]  = 0;
        stride[0] = FSPACE_DIM1;
        count[0]  = 1;
        block[0]  = FSPACE_DIM1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        /*
         * Create dataspace for the first dataset.
         */
        hsize_t mdim[datain.dimension()] = {FSPACE_DIM1};      /* Dimension size of the first dataset (in memory) */

        DataSpace mspace( (unsigned int) datain.dimension(), mdim );
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        dataset->write( datain.get_arr(), hdf5_type, mspace, fspace );
        delete dataset;
    }

    void writeStringdataset(string &str, string dsetname){
        const unsigned int FSPACE_DIM1 = str.size();
        hsize_t         str_dimsf[1] {FSPACE_DIM1};
        H5::DataSpace   dataspace(1, str_dimsf);

        DataSpace fspace( 1, str_dimsf );
        DataSpace mspace( 1, str_dimsf );

        hsize_t start[1];  // Start of hyperslab
        hsize_t stride[1]; // Stride of hyperslab
        hsize_t count[1];  // Block count
        hsize_t block[1];  // Block sizes
        start[0]  = 0;
        stride[0] = FSPACE_DIM1;
        count[0]  = 1;
        block[0]  = FSPACE_DIM1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        H5::StrType datatype(H5::PredType::C_S1, FSPACE_DIM1);
        H5::DataSet *str_dataset = new DataSet(file_->createDataSet(dsetname, datatype, dataspace));

        str_dataset->write(str, datatype, fspace, mspace);
        delete str_dataset;
    }

    template<typename T>
    void write2Ddataset(BootesArray<T> &datain, string name, PredType hdf5_type){
        if (datain.dimension() != 2){
            throw 1;
        }

        const unsigned int FSPACE_DIM1 = datain.shape()[0];
        const unsigned int FSPACE_DIM2 = datain.shape()[1];
        unsigned int FSPACE_RANK = 2;

        hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2}; // dim sizes of ds (on disk)
        DataSpace fspace( FSPACE_RANK, fdim );

        DataSet* dataset = new DataSet(file_->createDataSet(name, hdf5_type, fspace));
        /*
         * Select hyperslab for the dataset in the file, using 3x2 blocks,
         * (4,3) stride and (2,4) count starting at the position (0,1).
         */
        hsize_t start[2];  // Start of hyperslab
        hsize_t stride[2]; // Stride of hyperslab
        hsize_t count[2];  // Block count
        hsize_t block[2];  // Block sizes
        start[0]  = 0;           start[1]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2;
        block[0]  = FSPACE_DIM1; block[1]  = 1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        /*
         * Create dataspace for the first dataset.
         */
        hsize_t mdim[datain.dimension()] = {FSPACE_DIM1, FSPACE_DIM2};      /* Dimension size of the first dataset (in memory) */

        DataSpace mspace( (unsigned int) datain.dimension(), mdim );

        start[0]  = 0;           start[1]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2;
        block[0]  = FSPACE_DIM1; block[1]  = 1;
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        dataset->write( datain.get_arr(), hdf5_type, mspace, fspace );
        delete dataset;
    }

    template<typename T>
    void write3Ddataset(BootesArray<T> &datain, string name, PredType hdf5_type){
        if (datain.dimension() != 3){
            throw 1;
        }

        const unsigned int FSPACE_DIM1 = datain.shape()[0];
        const unsigned int FSPACE_DIM2 = datain.shape()[1];
        const unsigned int FSPACE_DIM3 = datain.shape()[2];
        unsigned int FSPACE_RANK = 3;

        hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3}; // dim sizes of ds (on disk)
        DataSpace fspace( FSPACE_RANK, fdim );

        DataSet* dataset = new DataSet(file_->createDataSet(name, hdf5_type, fspace));
        /*
         * Select hyperslab for the dataset in the file, using 3x2 blocks,
         * (4,3) stride and (2,4) count starting at the position (0,1).
         */
        hsize_t start[3];  // Start of hyperslab
        hsize_t stride[3]; // Stride of hyperslab
        hsize_t count[3];  // Block count
        hsize_t block[3];  // Block sizes
        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        /*
         * Create dataspace for the first dataset.
         */
        hsize_t mdim[datain.dimension()] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3};      /* Dimension size of the first dataset (in memory) */

        DataSpace mspace( (unsigned int) datain.dimension(), mdim );

        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        dataset->write( datain.get_arr(), hdf5_type, mspace, fspace );
        delete dataset;
    }

    template<typename T>
    void write4Ddataset(BootesArray<T> &datain, string name, PredType hdf5_type){
        if (datain.dimension() != 4){
            throw 1;
        }

        const unsigned int FSPACE_DIM1 = datain.shape()[0];
        const unsigned int FSPACE_DIM2 = datain.shape()[1];
        const unsigned int FSPACE_DIM3 = datain.shape()[2];
        const unsigned int FSPACE_DIM4 = datain.shape()[3];
        unsigned int FSPACE_RANK = 4;

        hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3, FSPACE_DIM4}; // dim sizes of ds (on disk)
        DataSpace fspace( FSPACE_RANK, fdim );

        DataSet* dataset = new DataSet(file_->createDataSet(name, hdf5_type, fspace));
        /*
         * Select hyperslab for the dataset in the file, using 3x2 blocks,
         * (4,3) stride and (2,4) count starting at the position (0,1).
         */
        hsize_t start[4];  // Start of hyperslab
        hsize_t stride[4]; // Stride of hyperslab
        hsize_t count[4];  // Block count
        hsize_t block[4];  // Block sizes
        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;           start[3]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;           stride[3] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3; count[3]  = FSPACE_DIM4;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;           block[3]  = 1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        /*
         * Create dataspace for the first dataset.
         */
        hsize_t mdim[datain.dimension()] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3, FSPACE_DIM4};      /* Dimension size of the first dataset (in memory) */

        DataSpace mspace( (unsigned int) datain.dimension(), mdim );

        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;           start[3]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;           stride[3] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3; count[3]  = FSPACE_DIM4;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;           block[3]  = 1;
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        dataset->write( datain.get_arr(), hdf5_type, mspace, fspace );
        delete dataset;
    }

    template<typename T>
    void write5Ddataset(BootesArray<T> &datain, string name, PredType hdf5_type){
        if (datain.dimension() != 5){
            throw 1;
        }

        const unsigned int FSPACE_DIM1 = datain.shape()[0];
        const unsigned int FSPACE_DIM2 = datain.shape()[1];
        const unsigned int FSPACE_DIM3 = datain.shape()[2];
        const unsigned int FSPACE_DIM4 = datain.shape()[3];
        const unsigned int FSPACE_DIM5 = datain.shape()[4];
        unsigned int FSPACE_RANK = 5;

        hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3, FSPACE_DIM4, FSPACE_DIM5}; // dim sizes of ds (on disk)
        DataSpace fspace( FSPACE_RANK, fdim );

        DataSet* dataset = new DataSet(file_->createDataSet(name, hdf5_type, fspace));
        /*
         * Select hyperslab for the dataset in the file, using 3x2 blocks,
         * (4,3) stride and (2,4) count starting at the position (0,1).
         */
        hsize_t start[5];  // Start of hyperslab
        hsize_t stride[5]; // Stride of hyperslab
        hsize_t count[5];  // Block count
        hsize_t block[5];  // Block sizes
        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;           start[3]  = 0;           start[4]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;           stride[3] = 1;           stride[4] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3; count[3]  = FSPACE_DIM4; count[4]  = FSPACE_DIM5;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;           block[3]  = 1;           block[4]  = 1;
        fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        /*
         * Create dataspace for the first dataset.
         */
        hsize_t mdim[datain.dimension()] = {FSPACE_DIM1, FSPACE_DIM2, FSPACE_DIM3, FSPACE_DIM4, FSPACE_DIM5};      /* Dimension size of the first dataset (in memory) */

        DataSpace mspace( (unsigned int) datain.dimension(), mdim );

        start[0]  = 0;           start[1]  = 0;           start[2]  = 0;           start[3]  = 0;           start[4]  = 0;
        stride[0] = FSPACE_DIM1; stride[1] = 1;           stride[2] = 1;           stride[3] = 1;           stride[4] = 1;
        count[0]  = 1;           count[1]  = FSPACE_DIM2; count[2]  = FSPACE_DIM3; count[3]  = FSPACE_DIM4; count[4]  = FSPACE_DIM5;
        block[0]  = FSPACE_DIM1; block[1]  = 1;           block[2]  = 1;           block[3]  = 1;           block[4]  = 1;
        mspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

        dataset->write( datain.get_arr(), hdf5_type, mspace, fspace );
        delete dataset;
    }

    void close(){
        file_->close();
    }

    // destructor
    ~Output(){
        delete file_;
    }
    private:
        string fn_;
        char mode_;
        H5File* file_;
};


class ReadOutput{
    public:
        string fn;
        H5File* file;
    ReadOutput(string fn_input){
        fn = fn_input;
        file = new H5File;
        *file = H5File( fn, H5F_ACC_RDONLY );

    }

    ~ReadOutput(){
        delete file;
    }
    template <typename T>
    void getAttribute(string att_name, T &value){
        Attribute att = file->openAttribute(att_name);
        DataType type = att.getDataType();
        att.read(type, &value);
    }

    template <typename T>
    T getAttribute(string att_name){
        T value;
        Attribute att = file->openAttribute(att_name);
        DataType type = att.getDataType();
        att.read(type, &value);
        return value;
    }

    string getString(string DataSetName){
        string data_out;
        std::string field_name(DataSetName);
        DataSet dataset = file->openDataSet(DataSetName);
        StrType datatype = dataset.getStrType();
        DataSpace dataspace = dataset.getSpace();

        std::string field_value;
        dataset.read(field_value, datatype, dataspace);

        return field_value;
    }

    template <typename T>
    void get4Ddata(string DataSetName, unsigned int hdf5Start[4], unsigned int hdf5Select[4], unsigned int outputshape[4], BootesArray<double> &data_out){

        DataSet dataset = file->openDataSet(DataSetName);
        H5T_class_t type_class = dataset.getTypeClass();

        //FloatType floatype = dataset.getFloatType();
        /*
        * Get order of datatype and print message if it's a little endian.
        */
        //H5std_string order_string;
        //H5T_order_t order = floatype.getOrder( order_string );
        // cout << order_string << endl;

        // size_t size = floatype.getSize();

        DataSpace dataspace = dataset.getSpace();
        /*
        * Get the number of dimensions in the dataspace.
        */

        int rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[4];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

        /*
        * Define hyperslab in the dataset; implicitly giving strike and
        * block NULL.
        */
        hsize_t      offset[4];      // hyperslab offset in the file
        hsize_t      count[4];       // size of the hyperslab in the file
        offset[0] = hdf5Start[0];               // starting location of hyperslab in 1st axis
        offset[1] = hdf5Start[1];               // starting location of hyperslab in 2nd axis
        offset[2] = hdf5Start[2];               // starting location of hyperslab in 3st axis
        offset[3] = hdf5Start[3];               // starting location of hyperslab in 4nd axis
        count[0]  = hdf5Select[0];               // number of locations to be selected in 1st axis
        count[1]  = hdf5Select[1];               // number of locations to be selected in 2nd axis
        count[2]  = hdf5Select[2];               // number of locations to be selected in 3st axis
        count[3]  = hdf5Select[3];               // number of locations to be selected in 4nd axis
        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset);

        /*
        * Define the memory dataspace.
        */
        hsize_t dimsm[4];
        dimsm[0] = outputshape[0];
        dimsm[1] = outputshape[1];
        dimsm[2] = outputshape[2];
        dimsm[3] = outputshape[3];
        unsigned int RANK_OUT = 4;
        DataSpace memspace( RANK_OUT, dimsm );
        /*
        * Define memory hyperslab.
        */
        hsize_t      offset_out[4];   // hyperslab offset in memory
        hsize_t      count_out[4];    // size of the hyperslab in memory
        offset_out[0] = 0;
        offset_out[1] = 0;
        offset_out[2] = 0;
        offset_out[3] = 0;
        count_out[0]  = outputshape[0];
        count_out[1]  = outputshape[1];
        count_out[2]  = outputshape[2];
        count_out[3]  = outputshape[3];
        memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );
        const int shape0 = outputshape[0];
        const int shape1 = outputshape[1];
        const int shape2 = outputshape[2];
        const int shape3 = outputshape[3];
        double *data_out_temp = new double[shape0*shape1*shape2*shape3];
        dataset.read( data_out_temp, PredType::NATIVE_DOUBLE, memspace, dataspace );
        /*
        for (int ss = 0; ss < outputshape[0]; ss ++){
            for (int kk = 0; kk < outputshape[1]; kk ++){
                for (int jj = 0; jj < outputshape[2]; jj ++){
                    for (int ii = 0; ii < outputshape[3]; ii ++){
                        data_out(ss, kk, jj, ii) = data_out_temp[((ss * shape2 + kk) * shape1 + jj) * shape0 + ii];
                    }
                }
            }
        }
        */
        for (int ii = 0; ii < shape0*shape1*shape2*shape3; ii ++){
            data_out(ii) = data_out_temp[ii];
        }
        delete []data_out_temp;
    }


    template <typename T>
    void get1Ddata(string DataSetName, unsigned int hdf5Start[1], unsigned int hdf5Select[1], unsigned int outputshape[1], BootesArray<double> &data_out){
        data_out.NewBootesArray(outputshape[0]);

        DataSet dataset = file->openDataSet(DataSetName);
        H5T_class_t type_class = dataset.getTypeClass();

        //FloatType floatype = dataset.getFloatType();
        /*
        * Get order of datatype and print message if it's a little endian.
        */
        //H5std_string order_string;
        //H5T_order_t order = floatype.getOrder( order_string );
        // cout << order_string << endl;

        // size_t size = floatype.getSize();

        DataSpace dataspace = dataset.getSpace();
        /*
        * Get the number of dimensions in the dataspace.
        */

        int rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[1];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

        /*
        * Define hyperslab in the dataset; implicitly giving strike and
        * block NULL.
        */
        hsize_t      offset[1];      // hyperslab offset in the file
        hsize_t      count[1];       // size of the hyperslab in the file
        offset[0] = hdf5Start[0];               // starting location of hyperslab in 1st axis
        count[0]  = hdf5Select[0];               // number of locations to be selected in 1st axis
        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset);

        /*
        * Define the memory dataspace.
        */
        hsize_t dimsm[1];
        dimsm[0] = outputshape[0];
        unsigned int RANK_OUT = 1;
        DataSpace memspace( RANK_OUT, dimsm );
        /*
        * Define memory hyperslab.
        */
        hsize_t      offset_out[1];   // hyperslab offset in memory
        hsize_t      count_out[1];    // size of the hyperslab in memory
        offset_out[0] = 0;
        count_out[0]  = outputshape[0];
        memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );
        const int shape0 = outputshape[0];
        double *data_out_temp = new double[shape0];
        dataset.read( data_out_temp, PredType::NATIVE_DOUBLE, memspace, dataspace );
        for (int ii = 0; ii < outputshape[0]; ii ++){
            data_out(ii) = data_out_temp[ii];
        }
        delete []data_out_temp;
    }

    void close(){
        file->close();
    }
};


#endif
