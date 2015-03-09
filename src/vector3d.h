#ifndef VECTOR3D
#define VECTOR3D

#include <cassert>
#include <cmath>
#include <string>

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

    const vector3d& operator = (const vector3d& vec) const {
        return vec;
    }

    double dot(const vector3d& vec){
        return (_data[0]*vec[0] + _data[1]*vec[1] + _data[2]*vec[2]);
    }

    vector3d cross(const vector3d& vec) const {
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

    double squared_norm() const{
        return (_data[0]*_data[0] + _data[1]*_data[1] + _data[2]*_data[2]);
    }

    double norm() const {
        return sqrt(squared_norm());
    }

    int size() const{
        return 3;
    }

    vector3d operator - (const vector3d& vec) const {
        vector3d result;
        result[0] = _data[0] - vec[0];
        result[1] = _data[1] - vec[1];
        result[2] = _data[2] - vec[2];
        return result;
    }

    vector3d operator + (const vector3d& vec) const {
        vector3d result;
        result[0] = _data[0] + vec[0];
        result[1] = _data[1] + vec[1];
        result[2] = _data[2] + vec[2];
        return result;
    }

    void print () {
        std::cout << _data[0] << "\t" << _data[1] << "\t" << _data[2] << std::endl;
    }

    vector3d operator / (const double& val) const {
        vector3d result;
        result[0] = _data[0] / val;
        result[1] = _data[1] / val;
        result[2] = _data[2] / val;
        return result;
    }

    vector3d operator * (const double& val) const {
        vector3d result;
        result[0] = _data[0] * val;
        result[1] = _data[1] * val;
        result[2] = _data[2] * val;
        return result;
    }

    void operator = (const double val) {
        _data[0] = val;
        _data[1] = val;
        _data[2] = val;
    }

    friend std::ostream& operator << (std::ostream& os, const vector3d& vec){
        os << std::scientific << vec[0] << "\t" << vec[1] << "\t" << vec[2];
        return os;
    }

    vector3d operator - () const {
        vector3d result;
        result[0] = -_data[0];
        result[1] = -_data[1];
        result[2] = -_data[2];
        return result;
    }

    operator vector3d      &()       { return static_cast<      vector3d&>(*this); }
    operator vector3d const&() const { return static_cast<const vector3d&>(*this); }

    double* begin(){
        return _data;
    }

    void operator += (const vector3d& vec) {
        _data[0] += vec[0];
        _data[1] += vec[1];
        _data[2] += vec[2];
    }

    void operator -= (const vector3d& vec) {
        _data[0] -= vec[0];
        _data[1] -= vec[1];
        _data[2] -= vec[2];
    }

};


#endif // VECTOR3D

