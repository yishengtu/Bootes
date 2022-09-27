#ifndef BootesArray_lastwork_HPP_
#define BootesArray_lastwork_HPP_
#include <string>
#include <memory>
#include <iostream>
#include <openacc.h>
//#include <cuda_runtime.h>

#ifdef __CUDACC__
    #define CUDA_CALLABLE_MEMBER __host__ __device__
#else
   #define CUDA_CALLABLE_MEMBER__host__ __device__
#endif

template<typename T>
class BootesArray {

    public:
    __host__ __device__ BootesArray(){
        arr_ = nullptr;
        shape_ = nullptr;
        allocated_ = false;
    }
    __host__ __device__ void updatehost(){ // update host copy of data
        //#pragma acc update self( arr_[0:arrsize_])
        //#pragma acc update self( shape_[0:dimension_] )
    }
    __host__ __device__ void updatedev(){ // update device copy of data
        //#pragma acc update device( arr_[0:arrsize_])
        //#pragma acc update device( shape_[0:dimension_])
    }
   __host__ __device__  __attribute__((nothrow)) void NewBootesArray(int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = 1;
            size3_ = 1;
            size4_ = 1;
            size5_ = 1;
            size6_ = 1;
            Allocate();
            shape_ = new int[1];
            shape_[0] = size1;
            dimension_ = 1;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])     // this[0, 1] is the class itself.
            // Attach: if array is created already, then attach pointer to the array.
            }
    __host__ __device__  __attribute__((nothrow)) void NewBootesArray(int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = 1;
            size4_ = 1;
            size5_ = 1;
            size6_ = 1;
            Allocate();
            shape_ = new int[2];
            shape_[1] = size1;
            shape_[0] = size2;
            dimension_ = 2;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])
            }
    __host__ __device__ __attribute__((nothrow)) void NewBootesArray(int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = 1;
            size5_ = 1;
            size6_ = 1;
            Allocate();
            shape_ = new int[3];
            shape_[2] = size1;
            shape_[1] = size2;
            shape_[0] = size3;
            dimension_ = 3;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])
            }
    __host__ __device__ __attribute__((nothrow)) void NewBootesArray(int size4, int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = size4;
            size5_ = 1;
            size6_ = 1;
            Allocate();
            shape_ = new int[4];
            shape_[3] = size1;
            shape_[2] = size2;
            shape_[1] = size3;
            shape_[0] = size4;
            dimension_ = 4;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])
            }
    __host__ __device__ __attribute__((nothrow)) void NewBootesArray(int size5, int size4, int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = size4;
            size5_ = size5;
            size6_ = 1;
            Allocate();
            shape_ = new int[5];
            shape_[4] = size1;
            shape_[3] = size2;
            shape_[2] = size3;
            shape_[1] = size4;
            shape_[0] = size5;
            dimension_ = 5;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])
            }
    __host__ __device__ __attribute__((nothrow)) void NewBootesArray(int size6, int size5, int size4, int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = size4;
            size5_ = size5;
            size6_ = size6;
            Allocate();
            shape_ = new int[6];
            shape_[5] = size1;
            shape_[4] = size2;
            shape_[3] = size3;
            shape_[2] = size4;
            shape_[1] = size5;
            shape_[0] = size6;
            dimension_ = 6;
            allocated_ = true;
            //#pragma acc enter data copyin(this[0:1]) create(arr_[0:arrsize_]) create(shape_[0:dimension_])
            }
    // destructor
    __host__ __device__ void clean(){
        delete[] arr_;
        delete[] shape_;
        allocated_ = false;
    }

    // destructor
    __host__ __device__ ~BootesArray(){
        clean();
    }

    // #pragma acc routine seq   // openACC // inline needed because in header files, if not inline multiple defs are created.
    __host__ __device__ inline T *data() {
	return arr_;
    }

    __host__ __device__ inline T &operator() (const int x1) {
        #ifdef DEBUG
            CHECKOK(x1);
        #endif
        return arr_[x1];
    }
    __host__ __device__ inline T &operator() (const int x2, const int x1) {
        #ifdef DEBUG
            CHECKOK(x1 + size1_*x2);
        #endif
        return arr_[x1 + size1_*x2];
    }
    __host__ __device__ inline T &operator() (const int x3, const int x2, const int x1) {
        #ifdef DEBUG
            CHECKOK(x1 + size1_*(x2 + size2_*x3));
        #endif
        return arr_[x1 + size1_*(x2 + size2_*x3)];
    }
    __host__ __device__ inline T &operator() (const int x4, const int x3, const int x2, const int x1) {
        #ifdef DEBUG
            CHECKOK(x1 + size1_*(x2 + size2_*(x3 + size3_ * x4)));
        #endif
        return arr_[x1 + size1_*(x2 + size2_*(x3 + size3_ * x4))];
    }
    __host__ __device__ inline T &operator() (const int x5, const int x4, const int x3, const int x2, const int x1) {
        #ifdef DEBUG
            CHECKOK(x1 + size1_*(x2 + size2_*(x3 + size3_ * (x4 + size4_ * x5))));
        #endif
        return arr_[x1 + size1_*(x2 + size2_*(x3 + size3_ * (x4 + size4_ * x5)))];
    }
    __host__ __device__ inline T &operator() (const int x6, const int x5, const int x4, const int x3, const int x2, const int x1) {
        #ifdef DEBUG
            CHECKOK(x1 + size1_*(x2 + size2_*(x3 + size3_ * (x4 + size4_ * (x5 + size5_ * x6)))));
        #endif
        return arr_[x1 + size1_*(x2 + size2_*(x3 + size3_ * (x4 + size4_ * (x5 + size5_ * x6))))];
    }

    //cloner
    __host__ __device__ BootesArray<T> &operator=(const BootesArray<T> &rhs){
        arr_ = new T[rhs.arrsize_];
        for (int ii = 0; ii < rhs.arrsize_; ii++){
            arr_[ii] = rhs.arr_[ii];
        }
        size1_ = rhs.size1_;
        size2_ = rhs.size2_;
        size3_ = rhs.size3_;
        size4_ = rhs.size4_;
        size5_ = rhs.size5_;
        size6_ = rhs.size6_;
        int shape_length = sizeof(rhs.shape_) / sizeof(rhs.shape_[0]);
        shape_ = new int[shape_length];
        for (int ii = 0; ii < shape_length; ii++){
            shape_[ii] = rhs.shape_[ii];
        }
        dimension_ = rhs.dimension_;
        arrsize_ = rhs.arrsize_;
        return *this;
    }

    __host__ __device__ inline int *shape(){
        return shape_;
    }

    __host__ __device__ inline int dimension(){
        return dimension_;
    }

    __host__ __device__ inline int arrsize(){
        return arrsize_;
    }

    __host__ __device__ T* get_arr(){
        return arr_;
    }

    __host__ __device__ inline bool checkallocated(){
        return allocated_;
    }

    __host__ __device__ inline void set_uniform(int x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (double) x;
        }
    }
    __host__ __device__ inline void set_uniform(float x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (double) x;
        }
    }
    __host__ __device__ inline void set_uniform(double x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (double) x;
        }
    }

    public:
        T *arr_;
        int size1_;
        int size2_;
        int size3_;
        int size4_;
        int size5_;
        int size6_;
        int *shape_;
        int dimension_;
        int arrsize_;
        bool allocated_;

    __host__ __device__ void Allocate(){
        arrsize_ = size1_*size2_*size3_*size4_*size5_*size6_;
        arr_ = new T[arrsize_];
    }


};



#endif
