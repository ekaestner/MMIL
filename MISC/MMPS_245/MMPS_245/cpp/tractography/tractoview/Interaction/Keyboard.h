#ifndef KEYBOARD_INCLUDED
#define KEYBOARD_INCLUDED

#include <vector>
#include <iostream.h>
#include <string>
#include "kqueue.h"

class Keyboard
{
public:
   enum Key
   {
      nul = 0, soh = 1, stx = 2, etx = 3, eot = 4, enq = 5, ack = 6, 
      BELL = 7,
      BACKSPACE = 8, 
      TAB = 9,
      NEWLINE = 10,
      VERTTAB = 11, 
      FORMFEED = 12,
      CR = 13,
      so = 14, si = 15, dle = 16,
      dc1 = 17, dc2 = 18, dc3 = 19, dc4 = 20, nak = 21, syn = 22, etb = 23, 
      can = 24, em = 25, sub = 26, 
      ESC = 27, 
      fs = 28, gs = 29, rs = 30, us = 31,
      SPACE = 32,
      BANG = 33,
      DOUBLEQUOTES = 34,
      HASH = 35,
      DOLLAR = 36,
      PERCENT = 37,
      AMPERSAND = 38,
      RIGHTQUOTE = '\'',
      LEFTPAREN = 40,
      RIGHTPAREN = 41,
      ASTERISK = '*',
      POS = '+', 
      COMMA = ',', 
      NEG = '-', 
      PERIOD = '.',
      SLASH = '/',
      ZERO = '0', ONE = '1', TWO = '2', THREE = '3', FOUR = '4', FIVE = '5', 
      SIX = '6', SEVEN = '7', EIGHT = '8', NINE = '9', 
      COLON = ':',
      SEMICOLON = ';',
      LESSTHAN = '<',
      EQUALS = '=',
      GREATERTHAN = '>',
      QUESTION = '?',
      AT = '@',
      A = 'A',B = 'B',C = 'C',D = 'D',E = 'E',F = 'F',G = 'G',H = 'H',
      I = 'I',J = 'J',K = 'K',L = 'L',M = 'M',N = 'N',O = 'O',P = 'P',Q = 'Q',
      R = 'R',S = 'S',T = 'T',U = 'U',V = 'V',W = 'W',X = 'X',Y = 'Y',Z = 'Z',
      LEFTBRACKET = '[', 
      BACKSLASH = '\\',
      RIGHTBRACKET = ']',
      CARET = '^',
      UNDERSCORE = '_',
      LEFTQUOTE = '`',
      a = 'a',b = 'b',c = 'c',d = 'd',e = 'e',f = 'f',g = 'g',h = 'h',
      i = 'i',j = 'j',k = 'k',l = 'l',m = 'm',n = 'n',o = 'o',p = 'p',q = 'q',
      r = 'r',s = 's',t = 't',u = 'u',v = 'v',w = 'w',x = 'x',y = 'y',z = 'z',
      CURLYLEFT = '{',
      PIPE = '|',
      CURLYRIGHT = '}',
      TILDE = '~',
      DEL = 127,
      /////////////////////////////////////////////////////////////////
      PAGEUP = 256, 
      PAGEDOWN = 257, 
      HOME = 258, 
      END = 259, 
      INSERT = 260,
      UPARROW = 261, 
      DOWNARROW = 262, 
      LEFTARROW = 263, 
      RIGHTARROW = 264,
      CAPS = 265, 
      LSHIFT = 266, 
      LCTRL = 267, 
      RSHIFT = 268, 
      RCTRL = 269, 
      F1 = 270, F2 = 271, F3 = 272, F4 = 273, F5 = 274, F6 = 275, 
      F7 = 276, F8 = 277, F9 = 278, F10 = 279, F11 = 280, F12 = 281,
      NumberOfKeys = 282 // make sure this number is bigger than the others.
   };
   enum BinaryState
   {
      ON, OFF
   };
   enum EdgeTriggerState
   {
	   DOWN, EDGEDOWN, EDGEUP, UP, FLOATING
   };
      
   Keyboard()
   {
      edges.resize( NumberOfKeys );
      binarys.resize( NumberOfKeys );
      for (int x = 0; x < NumberOfKeys; ++x)
      {
         edges[x] = UP;
         binarys[x] = OFF;
      }
   }
   
   const EdgeTriggerState& edgeState( char k ) const
   {
      int kk = (int)k;
      return edges[kk];
   }
   
   EdgeTriggerState& edgeState( char k )
   {
      int kk = (int)k;
      return edges[kk];
   }
   
   const EdgeTriggerState& edgeState( Key k ) const
   {
      assert( k >= 0 && k < edges.size() && "out of bounds" );
      return edges[k];
   }
   
   EdgeTriggerState& edgeState( Key k )
   {
      assert( k >= 0 && k < edges.size() && "out of bounds" );
      return edges[k];
   }
   
   const BinaryState& binaryState( char k ) const
   {
      int kk = (int)k;
      return binarys[kk];
   }
   
   BinaryState& binaryState( char k )
   {
      int kk = (int)k;
      return binarys[kk];
   }
   
   const BinaryState& binaryState( Key k ) const
   {
      assert( k >= 0 && k < binarys.size() && "out of bounds" );
      return binarys[k];
   }
   
   BinaryState& binaryState( Key k )
   {
      assert( k >= 0 && k < binarys.size() && "out of bounds" );
      return binarys[k];
   }
      
   void updateEdgeStates()
   {
      for (int x = 0; x < binarys.size(); ++x)
      {
         //find out if...
         switch( binarys[x] )
         {
         case ON:
            //button has been held down or it is just down.
            switch( edges[x] )
            {
            case FLOATING:
            case EDGEUP:
            case UP:       edges[x] = EDGEDOWN;  break;
            case EDGEDOWN: edges[x] = DOWN; break;
            }
         break;

         case OFF:
            //button has been up or it is just up.
            switch( edges[x] )
            {
            case FLOATING:
            case EDGEDOWN:
            case DOWN:    edges[x] = EDGEUP; break;
            case EDGEUP:  edges[x] = UP; break;
            }
         break;

         default: 
            //cout<<"broken on key: "<<x<<" of "<<binarys.size()<<" == "<<binarys[x]<<"\n"<<flush;
            assert(false);
         }
      }
   }

   
   kev::queue<Key> queue;
   std::vector<EdgeTriggerState> edges;
   std::vector<BinaryState> binarys;
};

#endif
