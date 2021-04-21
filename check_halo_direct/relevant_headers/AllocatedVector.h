// Author: Esteban Rangel
// erangel@anl.gov

#ifndef ALLOCATED_VECTOR_H
#define ALLOCATED_VECTOR_H
#include <stdexcept>
#include <algorithm>
#include <stdio.h>

template <class T>
class allocated_vector {
  private:
    size_t size_;
    size_t capacity_;
    T*  buff_;
  public:
    allocated_vector();
    allocated_vector(size_t, void*);
    ~allocated_vector();
    void set_buffer(void* buff) { buff_ = (T*)buff; }
    void set_capacity(size_t capacity) { capacity_ = capacity; }
    void push_back(T);
    void pop_back();
    size_t size();
    void resize(size_t);
    void resize(size_t, T);
    void shrink_to_fit();
    size_t capacity();
    size_t sizeof_item();
    T*   begin();
    T*   end();
    T*   cap_end();
    T&   at(size_t);
    T&   operator [](size_t i) { return buff_[i]; }
    allocated_vector& operator =(const allocated_vector&);
};

template <class T>
allocated_vector<T>::allocated_vector() {
  this->size_ = 0;
  this->capacity_ = 0;
  this->buff_ = NULL;
}

template <class T>
allocated_vector<T>::allocated_vector( size_t capacity, void* buff) {
  if (!buff)
    throw std::runtime_error("Cannot set buffer to NULL.");
  this->size_ = 0;
  this->capacity_ = capacity;
  this->buff_ = (T*)buff;
}

template <class T>
allocated_vector<T>::~allocated_vector() {
  ;
}

template <class T>
allocated_vector<T>& allocated_vector<T>::operator=( const allocated_vector& rhs ) {
  if ( capacity_ < rhs.size_ ) {
    char msg[256];
    sprintf( msg, "allocated_vector<T>::operator= Capacity exceeded [copy size=%lu capacity=%lu].", rhs.size_, capacity_ );
    throw std::runtime_error( (const char*)msg );
  }
  size_ = rhs.size_;
  std::copy( rhs.buff_, rhs.buff_+rhs.size_, buff_ );
  return *this;
}

template <class T>
size_t allocated_vector<T>::size() {
  return this->size_;
}

template <class T>
size_t allocated_vector<T>::capacity() {
  return this->capacity_;
}

template <class T>
void allocated_vector<T>::push_back( T val ) {
  if ( this->size_ + 1 > this->capacity_ )
    throw std::runtime_error("allocated_vector<T>::push_back( T val ) Capacity exceeded.");
  this->buff_[ this->size_++ ] = val;
}

template <class T>
void allocated_vector<T>::resize( size_t size ) {
  if ( size > this->capacity_ ) {
    char msg[256];
    sprintf( msg, "allocated_vector<T>::resize( size_t size ) Capacity exceeded [size=%lu capacity=%lu].", size, this->capacity_ );
    throw std::runtime_error( (const char*)msg );
  }
  this->size_ = size;
}

template <class T>
void allocated_vector<T>::resize( size_t size, T val ) {
  if ( size > this->capacity_ ) {
    char msg[256];
    sprintf( msg, "allocated_vector<T>::resize( size_t size, T val ) Capacity exceeded [size=%lu capacity=%lu].", size, this->capacity_ );
    throw std::runtime_error( (const char*)msg );
  }
  if ( this->size_ < size )
    std::fill( this->end(), this->end()+(size-this->size_), val );
  this->size_ = size;
}

template <class T>
void allocated_vector<T>::shrink_to_fit() {
  this->capacity_ = this->size_;
}

template <class T>
void allocated_vector<T>::pop_back() {
  if ( this->size > 0 )
    --this->size_;
}

template <class T>
T& allocated_vector<T>::at( size_t index ) {
  if ( index >= this->size_ )
    throw std::runtime_error("allocated_vector<T>::at( size_t index ) Out of range.");
  return this->buff_[ index ];
}

template <class T>
T* allocated_vector<T>::begin() {
  return this->buff_;
}

template <class T>
T* allocated_vector<T>::end() {
  return this->buff_+this->size_;
}

template <class T>
T* allocated_vector<T>::cap_end() {
  return this->buff_+this->capacity_;
}

template <class T>
size_t allocated_vector<T>::sizeof_item() {
  return sizeof(T);
}

#endif
