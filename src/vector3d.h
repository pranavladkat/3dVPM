#ifndef VECTOR3D
#define VECTOR3D

#include <cassert>
#include <cmath>

class vector3d{
private:
    double _data[3];
public:
    vector3d(){
        _data[0] = 0.0;
        _data[1] = 0.0;
        _data[2] = 0.0;
    }

    vector3d(double a,double b,double c){
        _data[0] = a;
        _data[1] = b;
        _data[2] = c;
    }

    vector3d(const vector3d& vec){
        _data[0] = vec[0];
        _data[1] = vec[1];
        _data[2] = vec[2];
    }

    double operator [] (int i) const{
        assert(i < 3);
        return _data[i];
    }

    double& operator [] (int i){
        assert(i < 3);
        return _data[i];
    }

    void operator = (const vector3d& vec){
        _data[0] = vec[0];
        _data[1] = vec[1];
        _data[2] = vec[2];
    }

    double dot(const vector3d& vec){
        return (_data[0]*vec[0] + _data[1]*vec[1] + _data[2]*vec[2]);
    }

    vector3d cross(const vector3d& vec){
        vector3d cross;
        cross[0] = _data[1]*vec[2] - _data[2]*vec[1];
        cross[1] = _data[2]*vec[0] - _data[0]*vec[2];
        cross[2] = _data[0]*vec[1] - _data[1]*vec[0];
        return cross;
    }

    void normalize(){
        double n = norm();
        _data[0] /= n;
        _data[1] /= n;
        _data[2] /= n;
    }

    double squared_norm(){
        return (_data[0]*_data[0] + _data[1]*_data[1] + _data[2]*_data[2]);
    }

    double norm(){
        return sqrt(squared_norm());
    }

    int size() const{
        return 3;
    }

};


#endif // VECTOR3D

