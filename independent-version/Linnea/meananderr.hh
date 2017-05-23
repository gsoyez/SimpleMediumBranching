#ifndef MEANANDERR_HH_INCLUDED
#define MEANANDERR_HH_INCLUDED

/// \class MeanAndErr
/// Get mean, error etc of a variable
///
///Each iteration, addEntry updates the values of mean etc

class MeanAndErr{
public:
    MeanAndErr();

    void addEntry(double myValue);
    double mean() const;
    double variance() const;
    double error() const;

private:
    unsigned long long _counter;
    double _value;
    double _valueSquared;
};

#endif // MEANANDERR_HH_INCLUDED
