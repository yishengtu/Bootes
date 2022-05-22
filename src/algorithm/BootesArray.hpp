#ifndef BootesArray_lastwork_HPP_
#define BootesArray_lastwork_HPP_
#include <string>
#include <memory>
#include <iostream>


template<typename T>
class BootesArray {

    public:
    BootesArray(){
        arr_ = nullptr;
        shape_ = nullptr;
        allocated_ = false;
    }
    //shape = {&size1_, &size2_, &size3_, &size4_, &size5_};
/*
    BootesArray(int size1) :
        size1(size1), size2(1),     size3(1),     size4(1),     size5(1)     { Allocate(); }
    BootesArray(int size1, int size2) :
        size1(size1), size2(size2), size3(1),     size4(1),     size5(1)     { Allocate(); }
    BootesArray(int size1, int size2, int size3) :
        size1(size1), size2(size2), size3(size3), size4(1),     size5(1)     { Allocate(); }
    BootesArray(int size1, int size2, int size3, int size4) :
        size1(size1), size2(size2), size3(size3), size4(size4), size5(1)     { Allocate(); }
    BootesArray(int size1, int size2, int size3, int size4, int size5) :
        size1(size1), size2(size2), size3(size3), size4(size4), size5(size5) { Allocate(); }
    */
    __attribute__((nothrow)) void NewBootesArray(int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = 1;
            size3_ = 1;
            size4_ = 1;
            size5_ = 1;
            Allocate();
            shape_ = new int[1];
            shape_[0] = size1;
            dimension_ = 1;
            allocated_ = true;
            }
    __attribute__((nothrow)) void NewBootesArray(int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = 1;
            size4_ = 1;
            size5_ = 1;
            Allocate();
            shape_ = new int[2];
            shape_[1] = size1;
            shape_[0] = size2;
            dimension_ = 2;
            allocated_ = true;
            }
    __attribute__((nothrow)) void NewBootesArray(int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = 1;
            size5_ = 1;
            Allocate();
            shape_ = new int[3];
            shape_[2] = size1;
            shape_[1] = size2;
            shape_[0] = size3;
            dimension_ = 3;
            allocated_ = true;
            }
    __attribute__((nothrow)) void NewBootesArray(int size4, int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = size4;
            size5_ = 1;
            Allocate();
            shape_ = new int[4];
            shape_[3] = size1;
            shape_[2] = size2;
            shape_[1] = size3;
            shape_[0] = size4;
            dimension_ = 4;
            allocated_ = true;
            }
    __attribute__((nothrow)) void NewBootesArray(int size5, int size4, int size3, int size2, int size1){
            if (allocated_){ clean(); }
            size1_ = size1;
            size2_ = size2;
            size3_ = size3;
            size4_ = size4;
            size5_ = size5;
            Allocate();
            shape_ = new int[5];
            shape_[4] = size1;
            shape_[3] = size2;
            shape_[2] = size3;
            shape_[1] = size4;
            shape_[0] = size5;
            dimension_ = 5;
            allocated_ = true;
            }
    // destructor
    void clean(){
        delete[] arr_;
        delete[] shape_;
        allocated_ = false;
    }

    // destructor
    ~BootesArray(){
        clean();
    }

    T &operator() (const int x1) {
        return arr_[x1];
    }
    T &operator() (const int x2, const int x1) {
        return arr_[x1 + size1_*x2];
    }
    T &operator() (const int x3, const int x2, const int x1) {
        return arr_[x1 + size1_*(x2 + size2_*x3)];
    }
    T &operator() (const int x4, const int x3, const int x2, const int x1) {
        return arr_[x1 + size1_*(x2 + size2_*(x3 + size3_ * x4))];
    }
    T &operator() (const int x5, const int x4, const int x3, const int x2, const int x1) {
        return arr_[x1 + size1_*(x2 + size2_*(x3 + size3_ * (x4 + size4_ * x5)))];
    }

    //cloner
    BootesArray<T> &operator=(const BootesArray<T> &rhs){
        arr_ = new T[rhs.arrsize_];
        for (int ii = 0; ii < rhs.arrsize_; ii++){
            arr_[ii] = rhs.arr_[ii];
        }
        size1_ = rhs.size1_;
        size2_ = rhs.size2_;
        size3_ = rhs.size3_;
        size4_ = rhs.size4_;
        size5_ = rhs.size5_;
        int shape_length = sizeof(rhs.shape_) / sizeof(rhs.shape_[0]);
        shape_ = new int[shape_length];
        for (int ii = 0; ii < shape_length; ii++){
            shape_[ii] = rhs.shape_[ii];
        }
        dimension_ = rhs.dimension_;
        arrsize_ = rhs.arrsize_;
        return *this;
    }

    int *shape(){
        return shape_;
    }

    int dimension(){
        return dimension_;
    }

    int arrsize(){
        return arrsize_;
    }

    T* get_arr(){
        return arr_;
    }

    void set_uniform(int x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (double) x;
        }
    }
    void set_uniform(float x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (float) x;
        }
    }
    void set_uniform(double x){
        for (int ii = 0; ii < arrsize_; ii ++){
            arr_[ii] = (double) x;
        }
    }

    private:
        T *arr_;
        int size1_;
        int size2_;
        int size3_;
        int size4_;
        int size5_;
        int *shape_;
        int dimension_;
        int arrsize_;
        bool allocated_;

    void Allocate(){
        arrsize_ = size1_*size2_*size3_*size4_*size5_;
        arr_ = new T[arrsize_];
    }
};


#endif


