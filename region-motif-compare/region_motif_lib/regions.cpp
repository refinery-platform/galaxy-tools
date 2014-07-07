#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

extern "C"  {
   typedef pair <int,int> Se_t;
   bool se_lt (const Se_t &l,const Se_t &r) { return l.first < r.first; }
   
   void merge_regions(int *regions, int*nregionsR,int *merge_sepR) {
      int nregs=nregionsR[0];
      if(nregs==0) return;
      int sep=merge_sepR[0];
      if(sep<1) sep=1;
      vector<Se_t> reg(nregs);
      for(int ireg=0;ireg<nregs;ireg++) {
	 reg[ireg]=make_pair(regions[ireg*2],regions[ireg*2+1]);
      }
      sort(reg.begin(),reg.end(),se_lt);
      int *reg_index = new int[nregs];
      for(int ireg=1;ireg<nregs;ireg++) reg_index[ireg]=-1;
      reg_index[0]=0;
      int last_ireg=0;
      int counter=1;
      for(int ireg=1;ireg<nregs;ireg++) {
	 if(reg[ireg].first<=reg[last_ireg].second+sep) {
	    if(reg[ireg].second>reg[last_ireg].second) reg[last_ireg].second=reg[ireg].second;
	 } else {
	    last_ireg=ireg;
	    reg_index[counter]=ireg;
	    counter++;
	 }
      }
      for(int ireg=0;ireg<counter;ireg++) {
	 regions[ireg*2] = reg[reg_index[ireg]].first;
	 regions[ireg*2+1] = reg[reg_index[ireg]].second;
      }
      nregionsR[0]=counter;
      delete [] reg_index;
   }
   void region_minus_region(int *regions, int*nregionsR,int *region2s, int*nregion2sR,int *updatedregions) {
      int sep=1;
      merge_regions(regions,nregionsR,&sep);
      merge_regions(region2s,nregion2sR,&sep);
      int nregs=nregionsR[0];
      int nreg2s=nregion2sR[0];
      for(int i=0;i<2*(nregs+nreg2s);i++) updatedregions[i]=-1;
      if(nregs==0) return;
      int ireg = 0;
      int iregout = 0;
      for(int ireg2=0; ireg2<nreg2s;ireg2++) {
	 if(ireg==nregs) break;
	 if(region2s[ireg2*2+1] < regions[2*ireg]) continue;
	 if(region2s[ireg2*2] > regions[2*ireg+1]) {
	    updatedregions[2*iregout] = regions[2*ireg];
	    updatedregions[2*iregout+1] = regions[2*ireg+1];
	    ireg++;
	    ireg2--;
	    iregout++;
	    continue;
	 }
	 int s = regions[ireg*2];
	 int e = regions[ireg*2+1];
	 int s2 = region2s[ireg2*2];
	 int e2 = region2s[ireg2*2+1];
	 if(s2<=s && e2>=e) {
	    ireg++;
	    ireg2--;
	 }
	 else if(s2<=s) {
	    regions[ireg*2] = e2+1;
	    continue;
	 } else if(e2>=e) {
	    updatedregions[2*iregout] = s;
	    updatedregions[2*iregout+1] = s2-1;
	    ireg2--;
	    iregout++;
	    ireg++;
	 } else {
	    updatedregions[2*iregout] = s;
	    updatedregions[2*iregout+1] = s2-1;
	    regions[ireg*2] = e2+1;
	    iregout++;
	    ireg2--;
	 }
      }
      while(ireg<nregs) {
	 updatedregions[2*iregout] = regions[2*ireg];
	 updatedregions[2*iregout+1] = regions[2*ireg+1];
	 ireg++;
	 iregout++;
      }
   }	    
   void intersection_of_regions(int *regions, int*nregionsR,int *region2s, int*nregion2sR,int *updatedregions) {
      int sep=1;
      merge_regions(regions,nregionsR,&sep);
      merge_regions(region2s,nregion2sR,&sep);
      int nregs=nregionsR[0];
      int nreg2s=nregion2sR[0];
      for(int i=0;i<2*(nregs+nreg2s);i++) updatedregions[i]=-1;
      if(nregs==0) return;
      if(nreg2s==0) return;
      int ireg2 = 0;
      int iregout = 0;
      for(int ireg=0; ireg<nregs;ireg++) {
	 if(ireg2==nreg2s) return;
	 if(regions[ireg*2+1] < region2s[2*ireg2]) continue;
	 if(regions[ireg*2] > region2s[2*ireg2+1]) {ireg2++; ireg--; continue;}
	 
	 int s = regions[ireg*2];
	 if(s<region2s[ireg2*2]) s = region2s[ireg2*2];
	 int e = regions[ireg*2+1];
	 if(e>region2s[ireg2*2+1]) {
	    e = region2s[ireg2*2+1];
	    ireg2++;
	    ireg--;
	 }
	 updatedregions[2*iregout] = s;
	 updatedregions[2*iregout+1] = e;
	 iregout++;
      }
   }	    
}
