/*
    Copyright 2013 Renwick James Hudspith

    This file (gramschmidt.h) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file gramschmidt.h
   @brief Perform a gram-schmidt reunitarisation of out gauge matrices 
   <a href="linkURL">http://en.wikipedia.org/wiki/Gram-Schmidt_process </a>
 */
#ifndef GLU_GRAMSCHMIDT_H
#define GLU_GRAMSCHMIDT_H

/**
   @fn void reunit2( double complex *__restrict U )
   @brief reunitarises U into itself
   @param U :: overwritten with a reunitarised version of itself
 */
void 
reunit2( double complex *__restrict U ) ;

/**
   @fn void latt_reunitU( struct site *__restrict lat )
   @brief reunitarise the whole lattice links
 */
void 
latt_reunitU( struct site *__restrict lat ) ;


#endif
