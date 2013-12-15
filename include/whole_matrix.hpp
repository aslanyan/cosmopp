#ifndef COSMO_PP_WHOLE_MATRIX_HPP
#define COSMO_PP_WHOLE_MATRIX_HPP

#include <vector>

/// Whole Matrix class.

/// This class represents a matrix in l-m space, which includes the off-diagonal elements.
class WholeMatrix
{
public:
    /// Constructor.
    
    /// Constructs a matrix with all elements equal to 0.
    /// \param lMin Minimum l in the matrix, must be non-negative.
    /// \param lMax Maximum l in the matrix, must be not less than lMin.
    WholeMatrix(int lMin, int lMax);
    
    /// Constructor.
    
    /// Constructs a matrix by reading from a file.
    /// \param fileName The name of the file.
    /// \param textFile Tells if it is a text file (true) or binary file (false).
    WholeMatrix(const char* fileName, bool textFile = false) { if(textFile) readFromTextFile(fileName); else readFromFile(fileName); }
    
    /// Copy-constructor.
    
    /// Constructs a whole matrix by identically copying another matrix.
    /// \param other The whole matrix to be copied.
    WholeMatrix(const WholeMatrix& other) : lMin_(other.lMin_), lMax_(other.lMax_), data_(other.data_) {}
    
    /// Allows reading of an element.
    /// \param l1 First index l.
    /// \param m1 First index m.
    /// \param l Second index l.
    /// \param m Second index m.
    /// \return The value of the element.
    double element(int l1, int m1, int l, int m) const;
    
    /// Allows reading and writing of an element.
    /// \param l1 First index l.
    /// \param m1 First index m.
    /// \param l Second index l.
    /// \param m Second index m.
    /// \return A reference to the element.
    double& element(int l1, int m1, int l, int m);
    
    /// Reads from a file in text format.
    /// \param fileName The name of the text file.
    void readFromTextFile(const char* fileName);
    
    /// Writes the matrix into a text file.
    /// \param fileName The name of the text file.
    void writeIntoTextFile(const char* fileName) const;
    
    /// Reads from a binary file. This is faster than reading from a text file. Binary files are also smaller in size.
    /// \param fileName The name of the binary file.
    void readFromFile(const char* fileName);
    
    /// Writes into a binary file. This is faster than writing into a text file. Binary files are also smaller in size.
    /// \param fileName The name of the binary file.
    void writeIntoFile(const char* fileName) const;
    
    /// Reads the minimum l.
    /// \return The minimum l.
    int getLMin() const { return lMin_; }
    
    /// Reads the maximum l.
    /// \return The maximum l.
    int getLMax() const { return lMax_; }
    
    //void rotate(double phi, double theta, double psi, bool firstIndexT = true, bool secondIndexT = true);
    
private:
    void initialize();
    bool checkIndices(int l, int m) const;
    
private:
    // first index is l', second index is l' + m', third index is l, fourth index is l + m
    typedef std::vector<std::vector<std::vector<std::vector<double> > > > DataType;
    
private:
    int lMin_;
    int lMax_;
    DataType data_;
};

#endif
