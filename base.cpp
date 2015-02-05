//
//  base.cpp
//  Project1_555
//
//  Created by Bing Hong Fu on 1/31/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#include <cstdio>
#include <stdexcept>
#include <cmath>
#include "base.h"

namespace gr {
    namespace trellis {
        
        bool
        dec2base(unsigned int num, int base, std::vector<int> &s)
        {
            int l = s.size();
            unsigned int n=num;
            for(int i=0;i<l;i++) {
                s[l-i-1] = n % base; //MSB first
                n /= base;
            }
            if(n!=0) {
                printf("Number %d requires more than %d digits.",num,l);
                return false;
            }
            else
                return true;
        }
        
        unsigned int
        base2dec(const std::vector<int> &s, int base)
        {
            int l = s.size();
            unsigned int num=0;
            for(int i=0;i<l;i++)
                num=num*base+s[i];
            return num;
        }
        
        bool
        dec2bases(unsigned int num, const std::vector<int> &bases, std::vector<int> &s)
        {
            int l = s.size();
            unsigned int n=num;
            for(int i=0;i<l;i++) {
                s[l-i-1] = n % bases[l-i-1];
                n /= bases[l-i-1];
            }
            if(n!=0) {
                printf("Number %d requires more than %d digits.",num,l);
                return false;
            }
            else
                return true;
        }
        
        unsigned int
        bases2dec(const std::vector<int> &s, const std::vector<int> &bases)
        {
            int l = s.size();
            unsigned int num=0;
            for(int i=0;i<l;i++)
                num = num * bases[i] + s[i];
            return num;
        }
        
    } /* namespace trellis */
} /* namespace gr */