/*

MIT License 
Copyright (c) 2023 Junyan Dai, Tobias Rubel, Yunheng Han, Erin Molloy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include<utility>
#include "bipartition.hpp"
#include<string>
#include<vector>
#include<algorithm>
#include<string.h>
#include<map>
#include<set>
#include<sstream>
#include "small_dollo_parsimony.h"
#include<cstdint>

//return tuple here
std::pair<unsigned int,unsigned int> score(std::string s1, std::string s2, std::string s3, int k, unsigned int d, std::vector<std::tuple<int,int,int,std::string>>& loc, std::unordered_map<std::string,int>& cnv_classes) {
  //std::get<0>(loc[i]) for chr
  //std::get<1>(loc[i]) for begin
  //std::get<2>(loc[i]) for end
  int n_linked = 0;
  int l_loss = 0;
  int l_local = 0;
  int r_loss = 0;
  int r_local = 0;
  int prev = -1;
  int chr = -1;
  int current_class = 1;
  
  for (int i = 0; i < k; i++) {
    if (s1[i] == '1' && s2[i] == '0') {
      
      l_loss++;
      //automatically checks i>0
      if(std::get<0>(loc[i]) == chr) { 
	if(prev == i-1 && (std::get<1>(loc[i]) - std::get<2>(loc[i-1])) < d) {
      	  l_local++;
	  n_linked++;
	  std::string cl = std::get<3>(loc[i]);
	  if(cnv_classes.find(cl) == cnv_classes.end()){
	  	cnv_classes[cl] = current_class;
	  }
	}
      }

      prev = i;
    }
    if (s1[i] == '1' && s3[i] == '0') {
     
      r_loss++;
      //automatically checks i>0
      if(std::get<0>(loc[i]) == chr) {
	if(prev == i-1 && (std::get<1>(loc[i]) - std::get<2>(loc[i-1])) < d){
	  r_local++;
	  n_linked++;
	  std::string cl = std::get<3>(loc[i]);
	  if(cnv_classes.find(cl) == cnv_classes.end()){
	  	cnv_classes[cl] = current_class;
	  }
	}
      }

      prev = i;
    }
    if(std::get<0>(loc[i]) != chr){
      current_class++;
      l_loss -= l_local;
      l_local = 0;
      r_loss -= r_local;
      r_local = 0;
      chr = std::get<0>(loc[i]);
    }

  }

  unsigned int loss = l_loss + r_loss;
  if(n_linked > 0){
  	return std::make_pair(loss,n_linked);
  }
  else{
  	return std::make_pair(loss,0);
  }
}

bool get_state_aux(Bipartition A, int i, uint8_t** C) {
  
  for (int j = 0; j < A.size(); j++) {
    if (C[i][j] == 1 && A.contain_index(j))
      return true;
  }

  return false;
  
}

bool all_missing_data(Bipartition A, int i, uint8_t** C) {
  int count = 0;
  for (int j = 0; j < A.size(); j++) {
    if (C[i][j] == 2 && A.contain_index(j)) {
      count++;
    }
  }

  return count == A.count();

}
// Input: A, B, k, C, S
// output: the state when B and A/B is the two children of A in the arbitrary opt tree using k characters.
std::string get_state(Bipartition A, Bipartition B, int k, uint8_t** C) {
  int V[k];
  boost::dynamic_bitset<> tbs(A.size());
  tbs.flip();
  Bipartition S(tbs); 
  memset(V, 0, sizeof(V));
  for (int i = 0; i < k; i++) {
    int flag = 0;

    if (all_missing_data(A, i, C)) {
      V[i] = 2;
      continue;
    }

    Bipartition D = B.complement(A);

    if (get_state_aux(B, i, C)) flag++;

    if (get_state_aux(D, i, C)) flag++;

    Bipartition E = A.complement(S);

    if (get_state_aux(E, i, C)) flag++;

    if (flag >= 2) V[i] = 1;
    
  }
  std::string R;
  R.reserve(k);
  for (int i = 0; i < k; i++)
    R.push_back(V[i] + '0');
  return R;

}
