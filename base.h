//
//  base.h
//  Project1_555
//
//  Created by Bing Hong Fu on 1/31/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#ifndef __Project1_555__base__
#define __Project1_555__base__

#include <vector>

namespace gr {
    namespace trellis {
        
        /*!
         * \brief  change base
         */
        
        bool dec2base(unsigned int num, int base, std::vector<int> &s);
        bool dec2bases(unsigned int num, const std::vector<int> &bases, std::vector<int> &s);
        unsigned int base2dec(const std::vector<int> &s, int base);
        unsigned int bases2dec(const std::vector<int> &s, const std::vector<int> &bases);
        
    } /* namespace trellis */
} /* namespace gr */
#endif /* defined(__Project1_555__base__) */
