#include "vec3.h"
#include <cmath>
#include <iostream>
using namespace std;
vec3::vec3()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

void vec3::print()
{
    // Will print matlab syntax vector. Output will be like: [2.09, 5.3, 9.1];
    cout << "[" << components[0] << ", " << components[1] << ", " << components[2] << "]" << endl;
}

vec3::vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

double vec3::dot(vec3 otherVector)
{
    return otherVector[0]*components[0] + otherVector[1]*components[1] + otherVector[2]*components[2];
}

vec3 vec3::cross(vec3 otherVector)
{
     return vec3(y()*otherVector.z()-z()*otherVector.y(), z()*otherVector.x()-x()*otherVector.z(),
                 x()*otherVector.y()-y()*otherVector.x());
}

double vec3::lengthSquared()
{
    return components[0]*components[0]
            + components[1]*components[1]
            + components[2]*components[2];
}

double vec3::length()
{
    return sqrt(lengthSquared());
}

vec3 &vec3::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    components[2] -= rhs;
    return *this;
}

vec3 &vec3::operator-=(vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

vec3 &vec3::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    components[2] += rhs;
    return *this;
}

vec3 &vec3::operator+=(vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

vec3 &vec3::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    components[2] *= rhs;
    return *this;
}

vec3 &vec3::operator*=(vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

vec3 &vec3::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    components[2] /= rhs;
    return *this;
}

vec3 &vec3::operator/=(vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

void vec3::print(string name)
{
    cout << name << " = ";
    print();
}

std::ostream &operator <<(std::ostream &os, const vec3 &myVector)
{
    os << "[" << myVector.x() << ", " << myVector.y()
       << ", " << myVector.z() << "];";
    return os;
}
