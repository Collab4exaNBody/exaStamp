#pragma once 
struct pass { template<typename ...T> pass(T...) {} };

template <class... T>
void unpackCopy(const char* p_buffer, size_t shift, size_t size, int nbElemPerMsg, T*... p_recv)
{
    assert( nbElemPerMsg >= size );
    size_t offset = 0;

    pass{(
        std::copy( (const T*)(&p_buffer[shift*sizeof(T) + offset]), (const T*)(&p_buffer[(shift+size)*sizeof(T) + offset]), p_recv) ,  
        offset += nbElemPerMsg*sizeof(T)
        ,1)...};

};


template <class... T>
void packCopy(const char* p_buffer, size_t shift, size_t begin, size_t end ,int nbElemPerMsg, T*... p_send)
{
    assert( nbElemPerMsg >= shift );
    size_t offset = 0;

    pass{(
        std::copy( p_send + begin, p_send + end, (T*) (&p_buffer[shift*sizeof(T) + offset]) ),
        offset += nbElemPerMsg*sizeof(T)
        ,1)...};

};

template <class U>
int sizeOfData(U* myData)
{
   return sizeof(U);  
};


template <class U, class... T>
int sizeOfData(U* myData, T*... data)
{
   return sizeof(U) + sizeOfData(data...) ; 
};
