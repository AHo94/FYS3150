#ifndef CELESTIALS_H
#define CELESTIALS_H


class Celestials
{
public:
    Celestials();
    void set_properties(double M, double R);
    double get_mass(void);
    double get_distance(void);
private:
    double mass;    // Mass of the celestial object
    double distance;    // Distance from the object to the Sun
};

#endif // CELESTIALS_H
