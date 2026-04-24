#ifndef INITIAL_DISTRIBUTIONS_HPP
#define INITIAL_DISTRIBUTIONS_HPP

#include <vector> 
#include <iostream> 
#include <cmath>

class distribution {
    public:
        virtual double operator() (double x, double y)=0;
        virtual void print()=0;
        virtual ~distribution() {}   
};

class w0_uniform : public distribution {
    public:
        w0_uniform();
        double operator() (double x, double y);
        void print();
};

class w0_zero : public distribution {
    public:
        w0_zero();
        double operator() (double x, double y);
        void print();
};

class w0_alfven : public distribution {
    public:
        w0_alfven(double kx_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
};

class w0_polarized_alfven : public distribution {
    public:
        w0_polarized_alfven(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double k_norm;
    double amp;
};

class w0_orszag_tang : public distribution {
    public:
        w0_orszag_tang(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};

class w0_orszag_tang_2 : public distribution {
    public:
        w0_orszag_tang_2(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};


class j0_uniform : public distribution {
    public:
        j0_uniform();

        double operator() (double x, double y);
        void print();
};

class j0_uniform_y : public distribution {
    public:
        j0_uniform_y();

        double operator() (double x, double y);
        void print();
};

class j0_current_sheet: public distribution {
    public:
        j0_current_sheet(double kx_j, double amp_j);

        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
};


class j0_alfven: public distribution {
    public:
        j0_alfven(double kx_j, double amp_j);

        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
};

class j0_polarized_alfven : public distribution {
    public:
        j0_polarized_alfven(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double k_norm;
    double amp;
};

class j0_orszag_tang : public distribution {
    public:
        j0_orszag_tang(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};

class j0_orszag_tang_2 : public distribution {
    public:
        j0_orszag_tang_2(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};


class combined_distribution : public distribution {
private:
    distribution* w_part;
    distribution* j_part;
    int sign_type;   // +1 -> q+ = w + j,  -1 -> q- = w - j

public:
    combined_distribution(distribution* w0, distribution* j0, int sign);
    double operator()(double x, double y);
    void print();
};

#endif
