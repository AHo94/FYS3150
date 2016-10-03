#ifndef VEC3_H
#define VEC3_H
#include <string>
#include <vector>

class vec3
{
public:
    vec3();
    vec3(const vec3&) = default;
    vec3(vec3&&) = default;
    vec3(double x, double y, double z);

    // Some vector properties
    double dot(vec3 otherVector);   // Dot product/scalar product
    vec3 cross(vec3 otherVector);   // Cross product

    double x() const { return components[0]; }
    double y() const { return components[1]; }
    double z() const { return components[2]; }
    void setX(double x) { components[0] = x; }
    void setY(double y) { components[1] = y; }
    void setZ(double z) { components[2] = z; }

    double &operator()(int index){return components[index];}    // Allows access to vector(0)
    double &operator[](int index){return components[index];}    // Allows access to vector[0]
    double lengthSquared();     // r^2 = x^2 + y^2 + z^2
    double length();            // r = sqrt(x^2 + y^2 + z^2), requires more computing time

    vec3 &operator-=(double rhs);  // Subtraction with a scalar
    vec3 &operator-=(vec3 rhs);    // Subtract componentwise
    vec3 &operator+=(double rhs);  // Addition with a scalar
    vec3 &operator+=(vec3 rhs);    // Addition componentwise
    vec3 &operator*=(double rhs);  // Multiplication with a scalar
    vec3 &operator*=(vec3 rhs);    // Multiplication componentwise
    vec3 &operator/=(double rhs);  // Dividing with a scalar
    vec3 &operator/=(vec3 rhs);    // Dividing componentwise


    // Implementing printing
    void print();
    void print(std::string name);
    friend std::ostream& operator << (std::ostream& os, const vec3& myVector);   // Allows cout << vector << endl;
private:
    double components[3];
};

inline vec3 operator+(vec3 lhs, double rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator+(double lhs, vec3 rhs) {
    rhs += lhs;
    return rhs;
}

inline vec3 operator+(vec3 lhs, vec3 rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator-(vec3 lhs, double rhs) {
    lhs -= rhs;
    return lhs;
}

inline vec3 operator-(double lhs, vec3 rhs) {
    rhs -= lhs;
    return rhs;
}

inline vec3 operator-(vec3 lhs, vec3 rhs) {
    lhs -= rhs;
    return lhs;
}


inline vec3 operator*(vec3 lhs, double rhs) {
    lhs *= rhs;
    return lhs;
}

inline vec3 operator*(double lhs, vec3 rhs) {
    rhs *= lhs;
    return rhs;
}

inline vec3 operator*(vec3 lhs, vec3 rhs) {
    lhs *= rhs;
    return lhs;
}


inline vec3 operator/(vec3 lhs, double rhs) {
    lhs /= rhs;
    return lhs;
}

inline vec3 operator/(double lhs, vec3 rhs) {
    rhs /= lhs;
    return rhs;
}

inline vec3 operator/(vec3 lhs, vec3 rhs) {
    lhs /= rhs;
    return lhs;
}

#endif // VEC3_H
