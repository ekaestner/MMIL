#ifndef CYCLIC_REDUNDANCY_CHECK_ISO3309
#define CYCLIC_REDUNDANCY_CHECK_ISO3309

/*
   CRCs are calculated using standard CRC methods with pre and post conditioning, 
   as defined by ISO 3309 [ISO-3309] or ITU-T V.42 [ITU-V42]. The CRC polynomial employed is 
   x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1

   The 32-bit CRC register is initialized to all 1's, and then the data from each byte 
   is processed from the least significant bit (1) to the most significant bit (128). 
   After all the data bytes are processed, the CRC register is inverted (its ones complement is taken).
   This value is transmitted (stored in the file) MSB first. For the purpose of separating into 
   bytes and ordering, the least significant bit of the 32-bit CRC is defined to be the 
   coefficient of the x^31 term. 

   Practical calculation of the CRC always employs a precalculated table to greatly accelerate 
   the computation. 
*/
#ifndef NULL
#define NULL 0L
#endif

class Crc32
{
public:
   Crc32();

   /* Return the CRC of the bytes buf[0..len-1]. */
   // NOTE: endianness matters, for example if your data is big endian, 
   // and you are on a little endian machine, you will need to reverse 
   // the ordering of the bytes before using this function.
   unsigned long crc(const unsigned char* const buf, int len);

protected:
   /* Update a running CRC with the bytes buf[0..len-1]--the CRC
      should be initialized to all 1's, and the transmitted value
      is the 1's complement of the final running CRC (see the
      crc() routine below)). */
   // NOTE: to begin a running crc, initialize crc to all 1s (0xffffffffL)
   unsigned long update_crc(unsigned long crc, const unsigned char* const buf, unsigned int len);


   // if DYNAMIC_CRC_TABLE was defined at compile time, then it checks whether 
   //      crc_table_empty is set and if not, it creates the lookup table.
   // otherwise it just returns the predefined lookup table.
   const unsigned long * get_crc_table();

   /* Make the table for a fast CRC. */
   void make_crc_table(void);
   
   /* Table of CRCs of all 8-bit messages. */
   static unsigned long crc_table[256];
   
   /* Flag: has the table been computed? Initially false. */
   static int crc_table_empty;
   
// alternate method...
protected:
   // alternate algorithm, for debugging...
   unsigned long update_crc_alt(unsigned long crc, const unsigned char* const buf, unsigned int len);
   // always uses a dynamic crc table (generated when needed rather than precomputed)
   const unsigned long* get_crc_table_alt();
   void make_crc_table_alt(void);
   static int crc_table_empty_alt;
   static unsigned long crc_table_alt[256];
};


class Crc32Adler
{
public:
   // largest prime smaller than 65536 
   static const unsigned long BASE;

   // NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1
   static const unsigned long NMAX;

   // for loop unrolling
   inline void DO1( unsigned long s1, unsigned long s2, const unsigned char* const buf, unsigned int i ) 
   {
      s1 += buf[i]; 
      s2 += s1;
   }
   inline void DO2(unsigned long s1, unsigned long s2, const unsigned char* const buf, unsigned int i)  
   {
      DO1(s1, s2, buf,i); 
      DO1(s1, s2, buf,i+1);
   }
   inline void DO4(unsigned long s1, unsigned long s2, const unsigned char* const buf, unsigned int i) 
   {
      DO2(s1, s2, buf,i); 
      DO2(s1, s2, buf,i+2);
   }
   inline void DO8(unsigned long s1, unsigned long s2, const unsigned char* const buf, unsigned int i)
   {
      DO4(s1, s2, buf,i); 
      DO4(s1, s2, buf,i+4);
   }
   inline void DO16(unsigned long s1, unsigned long s2, const unsigned char* const buf)   
   {
      DO8(s1, s2, buf,0);
      DO8(s1, s2, buf,8);
   }

   // get the CRC...
   // unsigned long c = crc( 0L, Z_NULL, 0 );
   //
   // while (read_buffer(buffer, length) != EOF) 
   // {
   //    c = crc( c, buffer, length);
   // }
   // if (c != original_crc) error();
   unsigned long crc( unsigned long adler, const unsigned char* const databuf, unsigned int length )
   {
       unsigned long s1 = adler & 0xffff;
       unsigned long s2 = (adler >> 16) & 0xffff;
       int k;
       const unsigned char* buf = databuf;

       if (buf == NULL) 
          return 1L;

       while (length > 0) {
           k = length < NMAX ? length : NMAX;
           length -= k;
           while (k >= 16) {
               DO16(s1, s2, buf);
	       buf += 16;
               k -= 16;
           }
           if (k != 0) do {
               s1 += *buf++;
	       s2 += s1;
           } while (--k);
           s1 %= BASE;
           s2 %= BASE;
       }
       return (s2 << 16) | s1;
   }
};

#endif
