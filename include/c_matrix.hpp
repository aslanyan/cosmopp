#ifndef COSMO_CPP_C_MATRIX_HPP
#define COSMO_CPP_C_MATRIX_HPP

#include <vector>
#include <string>

/// Covariance matrix in pixel space.

/// This class represents a covariance matrix in pixel space. By definition, the matrix is symmetric.
class CMatrix
{
public:
    /// Constructor.
    
    /// Initializes a matrix with all of the elements equal to 0.
    /// \param nPix The number of pixels.
    CMatrix(int nPix);
    
    /// Constructor.
    
    /// Initializes a matrix by reading from a binary file.
    /// \param fileName The name of the file.
    CMatrix(const char* fileName) { readFromFile(fileName); }
    
    /// Gives access to an element of the matrix that can be changed. Note that changing the element (i, j) automatically
    /// changes the element (j, i) too.
    /// \param i First index.
    /// \param j Second index.
    /// \return A reference to the element.
    double& element(int i, int j) { return matrix_[getIndex(i, j)]; }
    
    /// Allows reading of a given element.
    /// \param i First index.
    /// \param j Second index.
    /// \return The value of the element.
    double element(int i, int j) const {return matrix_[getIndex(i, j)]; }
    
    /// Reads the matrix from a binary file.
    /// \param fileName The name of the file.
    void readFromFile(const char* fileName);
    
    /// Reads the matrix from a text file.
    /// \param fileName The name of the file.
    void readFromTextFile(const char* fileName);
    
    /// Writes the matrix into a binary file.
    /// \param fileName The name of the file.
    void writeIntoFile(const char* fileName) const;
    
    /// Writes the matrix into a text file.
    /// \param fileName The name of the file.
    void writeIntoTextFile(const char* fileName) const;
    
    /// Returns the number of the pixels.
    /// \return The number of the pixels.
    int getNPix() const { return nPix_; }
    
    /// Allows reading of the comment.
    /// \return The comment string.
    const std::string& comment() const { return comment_; }
    
    /// Allows changing the comment.
    /// \return A reference to the comment string.
    std::string& comment() { return comment_; }
    
    /// Given a mask file, reduces the matrix to keep only the unmasked pixels.
    /// \param maskFileName The name of the file containing the mask (fits format).
    void maskMatrix(const char* maskFileName);
    
    /// Given a mask, reduces the matrix to keep only the unmasked pixels.
    /// \param goodPixels A vector containing the indices of the unmasked pixels.
    void maskMatrix(const std::vector<int>& goodPixels);
    
private:
    int getIndex(int i, int j) const;
    void initialize();
    
private:
    int nPix_;
    std::vector<double> matrix_;
    
    std::string comment_;
};

#endif
