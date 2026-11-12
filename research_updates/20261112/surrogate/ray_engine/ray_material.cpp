// Reciprocal, incoherent, opaque-surface geometric rays. NOT full-wave or a film model.
// Parameterized A0/A3 copied from the repaired core; all lengths are micrometres.
// A0 uses the original normalized cross-section curve remapped between endpoints.
// Remapped A0 members are not all claimed to be ideal CPC concentrators.
// Complex electric-field Jones vectors are transported in the global Cartesian frame.
// Convention: exp(-i omega t), passive N=n+i*k, p=s cross propagation_direction.
// Grey active faces use a specified absorptance and PEC-like residual reflection phase.
// Finite trace/cutoff residuals remain an explicit energy bound, never escaped energy.
// [ray_fem_check M3 extension] Surface M3: walls/land are Ti Fresnel; floor ("active") faces
// are a planar coherent single layer (film N1, thickness t) on semi-infinite Ti, evaluated at
// the local incidence cosine per local s/p (Airy). Floor loss is split into film (SiO2) and
// Ti-under-film using the power transmitted into the substrate. Lateral film edges and the
// Ti wall buried inside the film are NOT modelled (planar-stack floor only).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
constexpr double PI=3.1415926535897932384626433832795;
// EPS is conservative bounding-box padding only, not a ray travel distance.
constexpr double EPS=1e-8;
constexpr double INF=std::numeric_limits<double>::infinity();
using Z=std::complex<double>;
struct V {
  double x=0,y=0,z=0;
  V operator+(V b) const {return {x+b.x,y+b.y,z+b.z};}
  V operator-(V b) const {return {x-b.x,y-b.y,z-b.z};}
  V operator*(double a) const {return {a*x,a*y,a*z};}
  double operator[](int i) const {return i==0?x:(i==1?y:z);}
};
double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V cross(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm(V a){return std::sqrt(dot(a,a));}
V unit(V a){double n=norm(a);if(n==0)throw std::runtime_error("zero vector");return a*(1/n);}
V reflected(V d,V n){return unit(d-n*(2*dot(d,n)));}
struct CV {
  Z x=0,y=0,z=0;
  CV operator+(const CV& b) const {return {x+b.x,y+b.y,z+b.z};}
  CV operator*(double a) const {return {a*x,a*y,a*z};}
};
CV complex_vector(V a,Z b=1){return {a.x*b,a.y*b,a.z*b};}
Z dot(V a,const CV& b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double power(const CV& e){return std::norm(e.x)+std::norm(e.y)+std::norm(e.z);}
double coordinate_scale(V p){return std::max({1.0,std::abs(p.x),std::abs(p.y),std::abs(p.z)});}
double ray_t_min(V p){return 32*std::numeric_limits<double>::epsilon()*coordinate_scale(p);}

struct RNG {
  uint64_t s;
  uint64_t next(){uint64_t z=(s+=UINT64_C(0x9e3779b97f4a7c15));z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);return z^(z>>31);}
  double uniform(){return (next()>>11)*0x1.0p-53;}
};
struct Box {
  V lo={INF,INF,INF},hi={-INF,-INF,-INF};
  void add(V p){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
  void add(const Box& b){add(b.lo);add(b.hi);}
  bool hit(V p,V d,double cap) const {
    double a=0,b=cap;
    for(int k=0;k<3;++k){
      if(std::abs(d[k])<1e-15){if(p[k]<lo[k]-EPS||p[k]>hi[k]+EPS)return false;continue;}
      double c=(lo[k]-EPS-p[k])/d[k],e=(hi[k]+EPS-p[k])/d[k];if(c>e)std::swap(c,e);
      a=std::max(a,c);b=std::min(b,e);if(a>b)return false;
    }return b>ray_t_min(p);
  }
};
struct Tri {V a,e1,e2,n,center;Box box;bool active;double area;};
Tri triangle(V a,V b,V c,bool active){
  Tri t;t.a=a;t.e1=b-a;t.e2=c-a;t.n=unit(cross(t.e1,t.e2));t.center=(a+b+c)*(1.0/3);t.active=active;t.area=.5*norm(cross(t.e1,t.e2));t.box.add(a);t.box.add(b);t.box.add(c);return t;
}
double intersect(const Tri& t,V p,V d){
  V q=cross(d,t.e2);double det=dot(t.e1,q);if(std::abs(det)<1e-13)return INF;
  V v=p-t.a;double u=dot(v,q)/det;if(u< -1e-10||u>1+1e-10)return INF;
  V w=cross(v,t.e1);double b=dot(d,w)/det;if(b< -1e-10||u+b>1+1e-10)return INF;
  double dist=dot(t.e2,w)/det;return dist>ray_t_min(p)?dist:INF;
}
double plane_roundoff(const Tri& t,V p){
  double scale=std::max({coordinate_scale(p),coordinate_scale(t.a),coordinate_scale(t.a+t.e1),coordinate_scale(t.a+t.e2)});
  return 64*std::numeric_limits<double>::epsilon()*scale;
}
bool same_reflection_plane(const Tri& a,const Tri& b,V p){
  double tol=std::max(plane_roundoff(a,p),plane_roundoff(b,p));
  return std::abs(dot(a.n,b.n))>=1-64*std::numeric_limits<double>::epsilon()
      &&std::abs(dot(b.a-a.a,a.n))<=tol
      &&std::abs(dot(b.a+b.e1-a.a,a.n))<=tol
      &&std::abs(dot(b.a+b.e2-a.a,a.n))<=tol;
}
bool self_contact(const Tri& candidate,const Tri& previous,V p,V d,double t){
  if(!same_reflection_plane(candidate,previous,p))return false;
  double tol=std::max(plane_roundoff(candidate,p),plane_roundoff(previous,p));
  // A same-plane intersection caused only by the origin's roundoff error:
  // its distance can be amplified at grazing incidence, so use a plane bound.
  return std::abs(dot(p-candidate.a,candidate.n))<=tol
      && t*std::abs(dot(d,candidate.n))<=2*tol;
}
V intersection_point(const Tri& t,V p,V d){
  V q=cross(d,t.e2),v=p-t.a;double det=dot(t.e1,q);
  double u=dot(v,q)/det,b=dot(d,cross(v,t.e1))/det;
  // Reconstruct on the intersected triangle rather than advancing past it.
  return t.a+t.e1*u+t.e2*b;
}
struct Node {Box box;int left=-1,right=-1,first=0,count=0;};
struct Mesh {
  std::vector<Tri> tris;std::vector<int> order;std::vector<Node> nodes;
  void tri(V a,V b,V c,bool active=false){tris.push_back(triangle(a,b,c,active));}
  void quad(V a,V b,V c,V d,bool active=false){tri(a,b,c,active);tri(a,c,d,active);}
  int node(int first,int count){
    int id=static_cast<int>(nodes.size());nodes.emplace_back();Box box,centers;
    for(int i=first;i<first+count;++i){box.add(tris[order[i]].box);centers.add(tris[order[i]].center);}
    nodes[id].box=box;nodes[id].first=first;nodes[id].count=count;
    if(count>8){V ext=centers.hi-centers.lo;int axis=ext.y>ext.x?1:0;if(ext.z>ext[axis])axis=2;
      int mid=first+count/2;std::nth_element(order.begin()+first,order.begin()+mid,order.begin()+first+count,[&](int a,int b){return tris[a].center[axis]<tris[b].center[axis];});
      int l=node(first,count/2),r=node(mid,count-count/2);nodes[id].left=l;nodes[id].right=r;nodes[id].count=0;
    }return id;
  }
  void build(){order.resize(tris.size());std::iota(order.begin(),order.end(),0);nodes.reserve(2*tris.size());if(!tris.empty())node(0,static_cast<int>(tris.size()));}
  int closest(V p,V d,double& distance,int previous_face=-1) const {
    if(nodes.empty())return -1;
    int answer=-1;std::array<int,128> stack{};int n=1;stack[0]=0;
    while(n){const Node& b=nodes[stack[--n]];if(!b.box.hit(p,d,distance))continue;
      if(b.count){for(int j=b.first;j<b.first+b.count;++j){int i=order[j];double t=intersect(tris[i],p,d);if(t<distance){if(previous_face>=0&&self_contact(tris[i],tris[previous_face],p,d,t))continue;distance=t;answer=i;}}}
      else{if(n+2>static_cast<int>(stack.size()))throw std::runtime_error("BVH stack overflow");stack[n++]=b.left;stack[n++]=b.right;}
    }return answer;
  }
  double active_area() const {double a=0;for(const auto&t:tris)if(t.active)a+=t.area;return a;}
};
double cpc_f(double z){double s=std::sin(PI/12),c=std::cos(PI/12),K=20*(1+s);return(-s*(K+c*z)+std::sqrt(K*K+2*K*c*z))/(c*c)-10;}
std::vector<V> polygon(int sides,double radius,double z,double cx=0){
  std::vector<V> p;double phase=sides==6?PI/6:PI/4,circum=radius/std::cos(PI/sides);
  for(int k=0;k<sides;++k){double a=phase+2*PI*k/sides;p.push_back({cx+circum*std::cos(a),circum*std::sin(a),z});}return p;
}
struct Geometry {
  std::string name;Mesh mesh;double height=0,pitch=100,cell_area=10000,opening_area=0,bottom_apothem=10,upper_apothem=35;bool cavity=false,hex=false,flat=false;
  Geometry(const std::string& label,int segments,double requested_height=90,double bottom_width=20,double upper_width=2*cpc_f(90),double requested_pitch=80):name(label){
    if(name=="flat_black"||name=="flat_mirror"||name=="flat_patch"||name=="flat_ti"||name=="flat_active") {flat=true;opening_area=cell_area;return;}
    if(name=="A0"||name=="A1"||name=="A2"||name=="A3"){
      cavity=true;hex=true;height=requested_height;pitch=requested_pitch;bottom_apothem=bottom_width/2;upper_apothem=upper_width/2;
      if(!std::isfinite(height)||!std::isfinite(pitch)||!std::isfinite(bottom_apothem)||!std::isfinite(upper_apothem)||height<=0||pitch<=0||bottom_apothem<=0||bottom_apothem>=upper_apothem||upper_apothem>=pitch/2)throw std::runtime_error("invalid cavity dimensions or overlapping openings");
      cell_area=std::sqrt(3)*.5*pitch*pitch;double upper=upper_apothem;opening_area=2*std::sqrt(3)*upper*upper;
      double offset=name=="A1"?10:(name=="A2"?20:0);int m=name=="A3"?1:segments;
      auto previous=polygon(6,bottom_apothem,0,-offset);V fc={-offset,0,0};for(int k=0;k<6;++k)mesh.tri(fc,previous[k],previous[(k+1)%6],true);
      bool baseline=height==90&&bottom_width==20&&upper_width==2*cpc_f(90);
      for(int j=1;j<=m;++j){double z=height*j/m;
        double f=name=="A3"?bottom_apothem+(upper-bottom_apothem)*z/height:(baseline?cpc_f(z):bottom_apothem+(upper-bottom_apothem)*(cpc_f(90.0*j/m)-10)/(cpc_f(90)-10));
        double cx=-offset+offset*(f-bottom_apothem)/(upper-bottom_apothem);auto current=polygon(6,f,z,cx);
        for(int k=0;k<6;++k)mesh.quad(previous[k],previous[(k+1)%6],current[(k+1)%6],current[k]);
        previous=current;
      }
    }else if(name=="C"){
      cavity=true;height=40;opening_area=70*70;auto bottom=polygon(4,10,0,-15),top=polygon(4,35,height);V fc={-15,0,0};
      for(int k=0;k<4;++k){mesh.tri(fc,bottom[k],bottom[(k+1)%4],true);mesh.quad(bottom[k],bottom[(k+1)%4],top[(k+1)%4],top[k]);}
    }else if(name=="B"){
      height=40;auto bottom=polygon(4,35,0),top=polygon(4,10,height,-15);V tc={-15,0,height};
      for(int k=0;k<4;++k){mesh.tri(tc,top[k],top[(k+1)%4]);bool plus_x=(bottom[k].x>0&&bottom[(k+1)%4].x>0);mesh.quad(bottom[k],bottom[(k+1)%4],top[(k+1)%4],top[k],plus_x);}
    }else throw std::runtime_error("unknown structure");
    mesh.build();
    if(hex)validate_mesh();
  }
  void validate_mesh() const {
    using Vertex=std::tuple<long long,long long,long long>;
    std::map<std::pair<Vertex,Vertex>,int> edges;
    auto key=[](V p){return Vertex{std::llround(p.x*1e8),std::llround(p.y*1e8),std::llround(p.z*1e8)};};
    for(const auto& t:mesh.tris){
      if(!std::isfinite(t.area)||t.area<=1e-12||std::abs(norm(t.n)-1)>1e-12)throw std::runtime_error("degenerate mesh triangle");
      std::array<V,3> points={t.a,t.a+t.e1,t.a+t.e2};
      for(int j=0;j<3;++j){V p=points[j];if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||p.z< -1e-8||p.z>height+1e-8)throw std::runtime_error("mesh vertex outside height range");Vertex a=key(p),b=key(points[(j+1)%3]);if(b<a)std::swap(a,b);++edges[{a,b}];}
    }
    int opening_edges=0;
    for(const auto& item:edges){int count=item.second;if(count==1){const auto& edge=item.first;if(std::abs(std::get<2>(edge.first)*1e-8-height)>1e-7||std::abs(std::get<2>(edge.second)*1e-8-height)>1e-7)throw std::runtime_error("mesh crack away from opening");++opening_edges;}else if(count!=2)throw std::runtime_error("non-manifold mesh edge");}
    if(opening_edges!=6)throw std::runtime_error("mesh opening topology mismatch");
    double exact=2*std::sqrt(3)*bottom_apothem*bottom_apothem;
    if(std::abs(active_area()-exact)>1e-9*exact)throw std::runtime_error("mesh floor area mismatch");
  }
  double active_area() const {if(!flat)return mesh.active_area();if(name=="flat_black"||name=="flat_active")return cell_area;return name=="flat_patch"?2500:0;}
  V sample(RNG& r) const {double u=r.uniform()-.5,v=r.uniform()-.5;return hex?V{pitch*u+.5*pitch*v,std::sqrt(3)*.5*pitch*v,height}:V{pitch*u,pitch*v,height};}
  bool enter(V& p) const {
    if(!hex)return std::abs(p.x)<=35&&std::abs(p.y)<=35;
    double ap=upper_apothem;
    for(int i=-1;i<=1;++i){for(int j=-1;j<=1;++j){V center={pitch*i+.5*pitch*j,std::sqrt(3)*.5*pitch*j,0};V q=p-center;bool inside=true;
      for(int k=0;k<3;++k){double a=PI*k/3;if(std::abs(q.x*std::cos(a)+q.y*std::sin(a))>ap+1e-10){inside=false;break;}}
      if(inside){p=q;return true;}
    }}return false;
  }
};

struct Config {
  std::string structure="A0",surface="M1";
  double theta=0,phi=0,n=0,k=0,grey_alpha=.90,cutoff=1e-10,wavelength=0;
  double height=90,bottom_width=20,upper_width=2*cpc_f(90),pitch=80;
  // M3 floor film (planar single layer on semi-infinite Ti); lengths in micrometres.
  double film_n=0,film_k=0,film_thickness=0;
  uint64_t rays=10000,seed=1;int segments=128,max_bounces=128;
  bool self_test=false,fresnel_test=false,jones_probe=false,film_probe=false;
};
Config args(int argc,char**argv){
  Config c;
  for(int i=1;i<argc;++i){std::string a=argv[i];
    if(a=="--self-test"){c.self_test=true;continue;}if(a=="--fresnel"){c.fresnel_test=true;continue;}if(a=="--jones-probe"){c.jones_probe=true;continue;}if(a=="--film-probe"){c.film_probe=true;continue;}
    if(i+1>=argc)throw std::runtime_error("missing argument value");std::string v=argv[++i];
    if(a=="--structure")c.structure=v;else if(a=="--surface")c.surface=v;else if(a=="--theta")c.theta=std::stod(v);else if(a=="--phi")c.phi=std::stod(v);
    else if(a=="--n")c.n=std::stod(v);else if(a=="--k")c.k=std::stod(v);else if(a=="--grey-alpha")c.grey_alpha=std::stod(v);else if(a=="--cutoff")c.cutoff=std::stod(v);else if(a=="--wavelength-um")c.wavelength=std::stod(v);
    else if(a=="--height-um")c.height=std::stod(v);else if(a=="--bottom-width-um")c.bottom_width=std::stod(v);else if(a=="--upper-width-um")c.upper_width=std::stod(v);else if(a=="--pitch-um")c.pitch=std::stod(v);
    else if(a=="--film-n")c.film_n=std::stod(v);else if(a=="--film-k")c.film_k=std::stod(v);else if(a=="--film-thickness-um")c.film_thickness=std::stod(v);
    else if(a=="--rays")c.rays=std::stoull(v);else if(a=="--seed")c.seed=std::stoull(v);else if(a=="--segments")c.segments=std::stoi(v);else if(a=="--max-bounces")c.max_bounces=std::stoi(v);else throw std::runtime_error("unknown option: "+a);
  }
  if(c.self_test)return c;
  if(c.film_probe)c.surface="M3";
  if(c.surface!="M0"&&c.surface!="M1"&&c.surface!="M2"&&c.surface!="M3")throw std::runtime_error("surface must be M0, M1, M2, or M3");
  if(!std::isfinite(c.theta)||!std::isfinite(c.phi)||c.theta<0||c.theta>=90||c.rays==0||c.rays>1000000000ULL||c.segments<1||c.segments>4096||c.max_bounces<1||c.max_bounces>256)throw std::runtime_error("invalid geometry or sampling range");
  if(!std::isfinite(c.grey_alpha)||c.grey_alpha<0||c.grey_alpha>1||!std::isfinite(c.cutoff)||c.cutoff<0||c.cutoff>1e-10||!std::isfinite(c.wavelength)||c.wavelength<0)throw std::runtime_error("invalid alpha, cutoff, or wavelength");
  if(!std::isfinite(c.n)||!std::isfinite(c.k)||c.n<0||c.k<0)throw std::runtime_error("nonfinite or negative refractive index");
  bool need_material=c.surface!="M0"||c.structure=="flat_ti"||c.fresnel_test||c.jones_probe;
  if(need_material&&(!std::isfinite(c.n)||!std::isfinite(c.k)||c.n<=0||c.k<0||(c.k==0&&!c.fresnel_test)))throw std::runtime_error("opaque material requires finite n>0, k>0; --fresnel alone permits k=0");
  if(c.surface=="M3"&&(!std::isfinite(c.film_n)||!std::isfinite(c.film_k)||!std::isfinite(c.film_thickness)||c.film_n<=0||c.film_k<0||c.film_thickness<0||!(c.wavelength>0)))throw std::runtime_error("M3 requires --film-n>0, --film-k>=0, --film-thickness-um>=0 and --wavelength-um>0");
  return c;
}
// Ts,Tp: fraction of the incident local-s/p power transmitted into the semi-infinite substrate
// of a layered face; negative means a single opaque interface (all lost power in that face).
struct Coeff {Z rs,rp;double Ts=-1,Tp=-1;};
Coeff fresnel(double cosine,double n,double k){
  double c=std::max(0.0,std::min(1.0,cosine));Z m(n,k),m2=m*m,q=std::sqrt(m2-(1-c*c));
  // Passive normal wavevector: positive attenuation into the material.
  if(q.imag()<0||(q.imag()==0&&q.real()<0))q=-q;
  return {(c-q)/(c+q),(m2*c-q)/(m2*c+q)};
}
Z passive_q(Z eps,double s2){Z q=std::sqrt(eps-s2);if(q.imag()<0||(q.imag()==0&&q.real()<0))q=-q;return q;}
// Coherent single layer: air (eps0=1) | film eps1=N1^2, thickness d | semi-infinite substrate eps2=N2^2.
// Same sign conventions as fresnel(): s admittance Y=q, r_ij=(Y_i-Y_j)/(Y_i+Y_j);
// p admittance Y=eps/q, r_ij=(Y_j-Y_i)/(Y_j+Y_i) (H-field reflection; reproduces rp of fresnel()).
// r=(r01+r12 e^{2i beta})/(1+r01 r12 e^{2i beta}), beta=(2 pi/lambda) d q1, Im(q_j)>=0.
// Transmitted power into the substrate: s T=Re(q2)/q0 |t_E|^2, t_ij=2q_i/(q_i+q_j);
// p (H form, Z=q/eps) T=Re(q2/eps2)/q0 |t_H|^2, t_ij^H=2 Z_i/(Z_i+Z_j).
Coeff film_stack(double cosine,Z N1,double d_um,Z N2,double lambda_um){
  double c=std::max(0.0,std::min(1.0,cosine));
  if(c<=0){Coeff g{Z(-1),Z(-1)};g.Ts=0;g.Tp=0;return g;} // exactly grazing: total reflection
  double s2=1-c*c;Z e1=N1*N1,e2=N2*N2,q0(c,0),q1=passive_q(e1,s2),q2=passive_q(e2,s2);
  Z beta=(2*PI/lambda_um)*d_um*q1,ph=std::exp(Z(0,1)*beta),ph2=ph*ph;
  Z r01s=(q0-q1)/(q0+q1),r12s=(q1-q2)/(q1+q2),t01s=2.0*q0/(q0+q1),t12s=2.0*q1/(q1+q2);
  Z dens=1.0+r01s*r12s*ph2,rs=(r01s+r12s*ph2)/dens,ts=t01s*t12s*ph/dens;
  Z Y0=1.0/q0,Y1=e1/q1,Y2=e2/q2;
  Z r01p=(Y1-Y0)/(Y1+Y0),r12p=(Y2-Y1)/(Y2+Y1);
  Z Z0=q0,Z1=q1/e1,Z2=q2/e2,t01p=2.0*Z0/(Z0+Z1),t12p=2.0*Z1/(Z1+Z2);
  Z denp=1.0+r01p*r12p*ph2,rp=(r01p+r12p*ph2)/denp,tp=t01p*t12p*ph/denp;
  Coeff f{rs,rp};f.Ts=q2.real()/c*std::norm(ts);f.Tp=Z2.real()/c*std::norm(tp);
  return f;
}
Coeff film_coefficients(const Config& c,double cosine){return film_stack(cosine,Z(c.film_n,c.film_k),c.film_thickness,Z(c.n,c.k),c.wavelength);}
Coeff coefficients(const Config& c,bool active,double cosine,int force=0){
  // force 1=black, 2=perfect mirror, 3=Ti, used only by flat controls.
  if(force==1||(force==0&&c.surface=="M0"&&active))return {0,0};
  if(force==2||(force==0&&c.surface=="M0"))return {-1,1};
  if(force==0&&c.surface=="M2"&&active){double r=std::sqrt(1-c.grey_alpha);return {-r,r};}
  if(force==0&&c.surface=="M3"&&active)return film_coefficients(c,cosine);
  return fresnel(cosine,c.n,c.k);
}
struct Transform {
  V s,p_in,p_out,d_out;Coeff f;
  Transform(V d,V normal,const Config& c,bool active,int force=0){
    V normal_in=dot(d,normal)>0?normal*(-1):normal;
    d_out=reflected(d,normal_in);s=cross(d,normal_in);
    if(norm(s)<1e-12){V axis=std::abs(d.x)<=std::abs(d.y)&&std::abs(d.x)<=std::abs(d.z)?V{1,0,0}:(std::abs(d.y)<=std::abs(d.z)?V{0,1,0}:V{0,0,1});s=cross(d,axis);}
    s=unit(s);p_in=unit(cross(s,d));p_out=unit(cross(s,d_out));f=coefficients(c,active,-dot(d,normal_in),force);
  }
  CV apply(const CV& e) const {return complex_vector(s,f.rs*dot(s,e))+complex_vector(p_out,f.rp*dot(p_in,e));}
};
// FILM and TIFLOOR are sub-tallies of ACTIVE (layered M3 floor: film volume loss and
// Ti-under-film loss); LAND is the sub-tally of WALL deposited on the flat top land.
enum Metric {TOTAL=0,ACTIVE=1,WALL=2,ESCAPE=3,RESIDUAL=4,FILM=5,TIFLOOR=6,LAND=7,METRICS=8};
using Energy=std::array<double,METRICS>;
double max_split_error=0; // |(film+Ti_floor)-(before-remaining)| over all layered events
struct Result {
  std::array<Energy,2> energy{};std::array<CV,2> e;
  int bounces=0,reason=0;uint64_t wraps=0;
  // reason: 0 exhausted/escaped, 1 intensity cutoff, 2 bounce/event cap, 3 geometry unresolved.
  void finish(bool escape,int why){reason=why;for(int j=0;j<2;++j){energy[j][escape?ESCAPE:RESIDUAL]+=power(e[j]);energy[j][TOTAL]=energy[j][ACTIVE]+energy[j][WALL];}}
};
void interact(Result& r,V& d,V normal,bool active,const Config& c,int force=0,bool land=false){
  Transform tr(d,normal,c,active,force);
  for(int j=0;j<2;++j){double before=power(r.e[j]);CV after=tr.apply(r.e[j]);double remaining=power(after);
    if(!std::isfinite(remaining)||remaining>before+1e-12)throw std::runtime_error("non-passive or nonfinite reflection");
    if(remaining>before){after=after*std::sqrt(before/remaining);remaining=before;}
    double lost=before-remaining;
    r.energy[j][active?ACTIVE:WALL]+=lost;if(land)r.energy[j][LAND]+=lost;
    if(active&&tr.f.Ts>=0){
      // Planar stack: local s and p are orthogonal eigen-polarizations, so the power entering
      // the substrate is |E_s|^2 Ts + |E_p|^2 Tp; the film takes the remainder of the loss.
      double es=std::norm(dot(tr.s,r.e[j])),ep=std::norm(dot(tr.p_in,r.e[j]));
      double ti=es*tr.f.Ts+ep*tr.f.Tp;
      double film=es*(1-std::norm(tr.f.rs)-tr.f.Ts)+ep*(1-std::norm(tr.f.rp)-tr.f.Tp);
      max_split_error=std::max(max_split_error,std::abs(film+ti-lost));
      if(film< -1e-12||ti< -1e-12)throw std::runtime_error("negative layered absorption");
      r.energy[j][TIFLOOR]+=ti;r.energy[j][FILM]+=lost-ti; // exact bookkeeping: FILM+TIFLOOR==ACTIVE
    }
    r.e[j]=after;
  }
  d=tr.d_out;++r.bounces;
}
Result trace(const Geometry& g,V p,V d,const std::array<CV,2>& initial,const Config& c){
  Result r;r.e=initial;
  if(g.flat){bool active=g.name=="flat_black"||g.name=="flat_active"||(g.name=="flat_patch"&&std::abs(p.x)<=25&&std::abs(p.y)<=25);
    int force=g.name=="flat_black"?1:(g.name=="flat_mirror"?2:(g.name=="flat_ti"?3:0));
    interact(r,d,{0,0,1},active,c,force);r.finish(true,0);return r;
  }
  if(g.cavity){if(!g.enter(p)){interact(r,d,{0,0,1},false,c,0,true);r.finish(true,0);return r;}}
  else p.z+=4*ray_t_min(p); // Begin strictly above protrusion tops, outside solid.
  int previous_face=-1;
  uint64_t event_cap=static_cast<uint64_t>(c.max_bounces)*100+1000;
  for(uint64_t event=0;event<event_cap;++event){
    double ws=power(r.e[0]),wp=power(r.e[1]);
    if(ws==0&&wp==0){r.finish(false,0);return r;}
    if(ws<=c.cutoff&&wp<=c.cutoff){r.finish(false,1);return r;}
    if(r.bounces>=c.max_bounces){r.finish(false,2);return r;}
    double t_top=INF;if(d.z>1e-15)t_top=(g.height-p.z)/d.z;
    if(t_top<0&&d.z>0){r.finish(true,0);return r;}
    double nearest=t_top;int face=g.mesh.closest(p,d,nearest,previous_face);
    double t_floor=INF,t_x=INF,t_y=INF;
    if(!g.cavity){
      if(d.z< -1e-15)t_floor=-p.z/d.z;
      if(d.x>1e-15)t_x=(g.pitch*.5-p.x)/d.x;else if(d.x< -1e-15)t_x=(-g.pitch*.5-p.x)/d.x;
      if(d.y>1e-15)t_y=(g.pitch*.5-p.y)/d.y;else if(d.y< -1e-15)t_y=(-g.pitch*.5-p.y)/d.y;
      double seam=std::min(t_x,t_y);
      if(seam<nearest&&seam<t_floor){
        double tolerance=std::max(ray_t_min(p),32*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(seam)));
        if(seam< -tolerance){r.finish(false,3);return r;}p=p+d*std::max(0.0,seam);
        if(std::abs(t_x-seam)<=tolerance)p.x=-std::copysign(g.pitch*.5,d.x);
        if(std::abs(t_y-seam)<=tolerance)p.y=-std::copysign(g.pitch*.5,d.y);
        previous_face=-1; // A translated neighbouring cell has distinct faces.
        ++r.wraps;continue;
      }
      if(t_floor<nearest){if(t_floor< -ray_t_min(p)){r.finish(false,3);return r;}p=p+d*std::max(0.0,t_floor);p.z=0;interact(r,d,{0,0,1},false,c);previous_face=-1;continue;}
    }
    if(face>=0){const Tri&t=g.mesh.tris[face];p=intersection_point(t,p,d);interact(r,d,t.n,t.active,c);previous_face=face;continue;}
    if(std::isfinite(t_top)&&t_top>=0){r.finish(true,0);return r;}
    r.finish(false,3);return r;
  }
  r.finish(false,2);return r;
}
struct Stats {
  uint64_t count=0;std::array<long double,METRICS> mean{},m2{};
  void add(const Energy& x){++count;for(int k=0;k<METRICS;++k){long double delta=x[k]-mean[k];mean[k]+=delta/count;m2[k]+=delta*(x[k]-mean[k]);}}
  double se(int k) const {return count>1?std::sqrt(static_cast<double>(std::max(0.0L,m2[k])/(count-1)/count)):0;}
};
const char* names[METRICS]={"A_total","A_active","A_wall","escape","residual","A_film","A_Ti_floor","A_land"};
void print_vec(V a){std::cout<<'['<<a.x<<','<<a.y<<','<<a.z<<']';}
void print_complex(Z a){std::cout<<'['<<a.real()<<','<<a.imag()<<']';}
void print_vec(const CV& a){std::cout<<'[';print_complex(a.x);std::cout<<',';print_complex(a.y);std::cout<<',';print_complex(a.z);std::cout<<']';}
void print_stats(const Stats& s){
  std::cout<<'{';for(int k=0;k<METRICS;++k){if(k)std::cout<<',';std::cout<<'"'<<names[k]<<"\":"<<s.mean[k];}
  std::cout<<",\"absorption_upper_bound\":"<<s.mean[TOTAL]+s.mean[RESIDUAL]<<",\"closure\":"<<s.mean[TOTAL]+s.mean[ESCAPE]+s.mean[RESIDUAL]<<",\"se\":{";
  for(int k=0;k<METRICS;++k){if(k)std::cout<<',';std::cout<<'"'<<names[k]<<"\":"<<s.se(k);}std::cout<<"}}";
}
void fresnel_output(const Config& c){Coeff f=fresnel(std::cos(c.theta*PI/180),c.n,c.k);std::cout<<"{\"theta_deg\":"<<c.theta<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"rs\":";print_complex(f.rs);std::cout<<",\"rp\":";print_complex(f.rp);std::cout<<",\"Rs\":"<<std::norm(f.rs)<<",\"Rp\":"<<std::norm(f.rp)<<",\"p_basis\":\"s_cross_d\"}\n";}
void film_output(const Config& c){
  double cs=std::cos(c.theta*PI/180);Coeff f=film_coefficients(c,cs),b=fresnel(cs,c.n,c.k);
  std::cout<<"{\"theta_deg\":"<<c.theta<<",\"wavelength_um\":"<<c.wavelength<<",\"film_n\":"<<c.film_n<<",\"film_k\":"<<c.film_k<<",\"film_thickness_um\":"<<c.film_thickness<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"rs\":";print_complex(f.rs);std::cout<<",\"rp\":";print_complex(f.rp);
  std::cout<<",\"Rs\":"<<std::norm(f.rs)<<",\"Rp\":"<<std::norm(f.rp)<<",\"Ts_substrate\":"<<f.Ts<<",\"Tp_substrate\":"<<f.Tp<<",\"As_film\":"<<1-std::norm(f.rs)-f.Ts<<",\"Ap_film\":"<<1-std::norm(f.rp)-f.Tp;
  std::cout<<",\"bare_rs\":";print_complex(b.rs);std::cout<<",\"bare_rp\":";print_complex(b.rp);std::cout<<",\"p_basis\":\"s_cross_d\",\"p_reflection\":\"H_field_equal_to_fresnel_rp_convention\"}\n";
}
void jones_output(const Config& c){
  V d=unit({.3,-.4,-.8});V s=unit(cross(d,{0,0,1})),p=unit(cross(s,d));
  CV e=complex_vector(s,std::sqrt(.6))+complex_vector(p,Z(.3,.4)*std::sqrt(1.6));
  std::array<V,2> normals={unit({.1,.2,1}),unit({.7,-.3,-1})};
  std::cout<<"{\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"initial_direction\":";print_vec(d);std::cout<<",\"initial_E\":";print_vec(e);std::cout<<",\"normals\":[";print_vec(normals[0]);std::cout<<',';print_vec(normals[1]);std::cout<<"],\"steps\":[";
  for(int j=0;j<2;++j){if(j)std::cout<<',';double before=power(e);Transform tr(d,normals[j],c,false,3);e=tr.apply(e);d=tr.d_out;std::cout<<"{\"direction\":";print_vec(d);std::cout<<",\"E\":";print_vec(e);std::cout<<",\"power_before\":"<<before<<",\"power_after\":"<<power(e)<<'}';}std::cout<<"]}\n";
}
std::array<CV,2> initial_fields(double phi,V d){V s={-std::sin(phi),std::cos(phi),0};return {complex_vector(s),complex_vector(unit(cross(s,d)))};}
void self_test(){
  Config c;c.n=2;c.k=3;c.surface="M1";V d={0,0,-1};auto initial=initial_fields(0,d);Geometry flat("flat_ti",1);
  auto out=trace(flat,{0,0,0},d,initial,c);double exact=1-std::norm((Z(1)-Z(2,3))/(Z(1)+Z(2,3)));
  if(std::abs(out.energy[0][TOTAL]-exact)>1e-13||std::abs(out.energy[1][TOTAL]-exact)>1e-13)throw std::runtime_error("normal Fresnel failed");
  Transform normal(d,{0,0,1},c,false);Coeff f=fresnel(1,2,3);CV ex=normal.apply(complex_vector({1,0,0}));if(std::abs(ex.x-f.rs)>1e-13)throw std::runtime_error("normal global polarization failed");
  Geometry a("A0",128),b("B",1),cavity("C",1);
  if(std::abs(a.active_area()-200*std::sqrt(3))>1e-8||std::abs(b.active_area()-1800*std::sqrt(2))>1e-8||std::abs(cavity.active_area()-400)>1e-8)throw std::runtime_error("geometry areas failed");
  RNG rng{7};for(const char* surface:{"M0","M1","M2"}){c.surface=surface;for(const Geometry* g:{&a,&b,&cavity}){for(int j=0;j<100;++j){auto r=trace(*g,g->sample(rng),d,initial,c);for(const auto& en:r.energy)if(std::abs(en[TOTAL]+en[ESCAPE]+en[RESIDUAL]-1)>1e-11)throw std::runtime_error("energy closure failed");}}}
  c.surface="M2";Geometry active("flat_active",1);auto grey=trace(active,{0,0,0},d,initial,c);if(std::abs(grey.energy[0][ACTIVE]-.9)>1e-13)throw std::runtime_error("grey active failed");
  // --- M3 layered-floor checks (ray_fem_check extension) ---
  // (a) zero thickness reproduces the bare Fresnel coefficients and puts all loss in the substrate.
  for(double th:{0.0,15.0,30.0,45.0,60.0,75.0,85.0,89.0}){double cs=std::cos(th*PI/180);
    for(Z N1:{Z(2.2447212969283283,0.06982232081911283),Z(0.5947997319778178,1.7674602002464557),Z(1,0)}){
      Coeff f=film_stack(cs,N1,0.0,Z(4.308219763636363,20.63042690909091),10.5),b=fresnel(cs,4.308219763636363,20.63042690909091);
      if(std::abs(f.rs-b.rs)>1e-12||std::abs(f.rp-b.rp)>1e-12||std::abs(f.Ts-(1-std::norm(b.rs)))>1e-12||std::abs(f.Tp-(1-std::norm(b.rp)))>1e-12)throw std::runtime_error("zero-thickness film does not reproduce bare Fresnel");}}
  // (b) passivity of the stack and of the film share over angle and thickness.
  for(double t:{0.1,1.0,2.0,3.0,6.0})for(double th=0;th<90;th+=2.5){Coeff f=film_stack(std::cos(th*PI/180),Z(2.2447212969283283,0.06982232081911283),t,Z(4.308219763636363,20.63042690909091),10.5);
    if(f.Ts< -1e-14||f.Tp< -1e-14||1-std::norm(f.rs)-f.Ts< -1e-12||1-std::norm(f.rp)-f.Tp< -1e-12)throw std::runtime_error("film stack not passive");}
  // (c) a lossless film (k=0) on Ti: all floor loss must be in the substrate.
  {Coeff f=film_stack(std::cos(40*PI/180),Z(1.5,0),2.0,Z(4.308219763636363,20.63042690909091),10.5);if(std::abs(1-std::norm(f.rs)-f.Ts)>1e-12||std::abs(1-std::norm(f.rp)-f.Tp)>1e-12)throw std::runtime_error("lossless film absorbs");}
  // (d) flat M3 face at normal incidence: tallies equal the analytic stack values.
  c.surface="M3";c.n=4.308219763636363;c.k=20.63042690909091;c.film_n=2.2447212969283283;c.film_k=0.06982232081911283;c.film_thickness=2;c.wavelength=10.5;
  {auto out3=trace(active,{0,0,0},d,initial,c);Coeff f=film_stack(1,Z(c.film_n,c.film_k),2,Z(c.n,c.k),10.5);
    for(int j=0;j<2;++j){double A=1-(j?std::norm(f.rp):std::norm(f.rs)),T=j?f.Tp:f.Ts;
      if(std::abs(out3.energy[j][ACTIVE]-A)>1e-13||std::abs(out3.energy[j][TIFLOOR]-T)>1e-13||std::abs(out3.energy[j][FILM]-(A-T))>1e-13)throw std::runtime_error("flat M3 tallies failed");}}
  // (e) M3 energy closure and sub-tally consistency in an A3 cavity at oblique incidence.
  {Geometry a3("A3",1,77,16.4,69.5455842631299,80);V dd=unit({-.35,-.2,-.91});auto init=initial_fields(std::atan2(.2,.35),dd);
    for(int j=0;j<300;++j){auto r=trace(a3,a3.sample(rng),dd,init,c);for(const auto& en:r.energy){
      if(std::abs(en[TOTAL]+en[ESCAPE]+en[RESIDUAL]-1)>1e-11)throw std::runtime_error("M3 energy closure failed");
      if(std::abs(en[FILM]+en[TIFLOOR]-en[ACTIVE])>1e-13||en[LAND]>en[WALL]+1e-15||en[FILM]< -1e-13||en[TIFLOOR]< -1e-13)throw std::runtime_error("M3 sub-tally inconsistency");}}}
  std::cout<<"{\"self_test\":\"PASS\",\"checks\":[\"normal_Fresnel\",\"normal_global_Jones_sign\",\"active_geometry_areas\",\"M0_M1_M2_energy_closure_900_rays\",\"grey_absorptance\",\"M3_zero_thickness_equals_bare_Fresnel_8_angles_3_films\",\"M3_stack_passivity_5_thicknesses_36_angles\",\"M3_lossless_film_no_film_loss\",\"M3_flat_normal_tallies\",\"M3_A3_energy_closure_and_subtallies_300_rays\"],\"max_split_error\":"<<max_split_error<<"}\n";
}
} // namespace
int main(int argc,char**argv){try{
  Config c=args(argc,argv);std::cout<<std::setprecision(17);
  if(c.self_test){self_test();return 0;}if(c.fresnel_test){fresnel_output(c);return 0;}if(c.jones_probe){jones_output(c);return 0;}if(c.film_probe){film_output(c);return 0;}
  auto began=std::chrono::steady_clock::now();Geometry g(c.structure,c.segments,c.height,c.bottom_width,c.upper_width,c.pitch);RNG rng{c.seed};double th=c.theta*PI/180,ph=c.phi*PI/180;
  V d={-std::sin(th)*std::cos(ph),-std::sin(th)*std::sin(ph),-std::cos(th)};auto initial=initial_fields(ph,d);
  std::array<Stats,3> stats;uint64_t bounces=0,wraps=0,cutoffs=0,caps=0,lost=0;int max_seen=0;double max_closure_error=0;
  for(uint64_t i=0;i<c.rays;++i){Result r=trace(g,g.sample(rng),d,initial,c);Energy avg{};
    for(int j=0;j<2;++j){stats[j].add(r.energy[j]);for(int k=0;k<METRICS;++k)avg[k]+=.5*r.energy[j][k];max_closure_error=std::max(max_closure_error,std::abs(r.energy[j][TOTAL]+r.energy[j][ESCAPE]+r.energy[j][RESIDUAL]-1));}
    stats[2].add(avg);bounces+=r.bounces;wraps+=r.wraps;max_seen=std::max(max_seen,r.bounces);cutoffs+=r.reason==1;caps+=r.reason==2;lost+=r.reason==3;
  }
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  std::cout<<"{\"schema\":1,\"model\":\"OPAQUE_POLARIZED_GEOMETRIC_RAYS_NOT_FULL_WAVE\",\"surface\":\""<<c.surface<<"\",\"structure\":\""<<c.structure<<"\",\"theta_deg\":"<<c.theta<<",\"phi_deg\":"<<c.phi<<",\"wavelength_um\":"<<c.wavelength<<",\"n\":"<<c.n<<",\"k\":"<<c.k<<",\"grey_alpha\":"<<c.grey_alpha<<",\"cutoff\":"<<c.cutoff<<",\"rays\":"<<c.rays<<",\"seed\":"<<c.seed<<",\"curve_segments\":"<<c.segments<<",\"max_bounces\":"<<c.max_bounces<<",\"channels\":{\"s\":";
  print_stats(stats[0]);std::cout<<",\"p\":";print_stats(stats[1]);std::cout<<",\"unpolarized\":";print_stats(stats[2]);
  std::cout<<"},\"ray_origin_policy\":\"exact_surface_barycentric_scale_roundoff_self_contact_v2\",\"max_ray_closure_error\":"<<max_closure_error<<",\"cutoff_terminated_rays\":"<<cutoffs<<",\"cap_terminated_rays\":"<<caps<<",\"geometry_unresolved_rays\":"<<lost<<",\"mean_bounces\":"<<double(bounces)/c.rays<<",\"max_observed_bounces\":"<<max_seen<<",\"periodic_wraps\":"<<wraps<<",\"wall_seconds\":"<<seconds<<",\"unit_cell_area_um2\":"<<g.cell_area<<",\"active_area_um2\":"<<g.active_area()<<",\"opening_area_um2\":"<<g.opening_area<<",\"triangles\":"<<g.mesh.tris.size()<<",\"height_um\":"<<c.height<<",\"bottom_width_um\":"<<c.bottom_width<<",\"upper_width_um\":"<<c.upper_width<<",\"pitch_um\":"<<c.pitch<<",\"mesh_validation\":\"PASS_manifold_except_six_opening_edges_positive_area_finite_bounds\",\"curve_family\":\"original_normalized_section_remapped_between_endpoints_not_all_ideal_CPC\",\"lattice\":\""<<(g.hex?"triangular_parameterized_primitive":"square_100um")<<"\",\"geometry_units\":\"um\",\"polarization\":\"paired_s_p_complex_global_Jones\",\"grey_residual_phase\":\"PEC_like_rs_negative_rp_positive\",\"uncertainty\":\"MC_SE_for_spatial_sampling_only_residual_is_separate_bound\",\"film_n\":"<<c.film_n<<",\"film_k\":"<<c.film_k<<",\"film_thickness_um\":"<<c.film_thickness<<",\"max_layered_split_error\":"<<max_split_error<<",\"limitations\":\""<<(c.surface=="M3"?"M3_planar_coherent_single_layer_floor_on_semi_infinite_Ti_local_cosine_no_lateral_film_edges_no_buried_wall_no_diffraction_or_cavity_interference":"opaque_local_bulk_Fresnel_no_diffraction_interference_roughness_film_or_temperature_model")<<"\"}\n";
  return 0;
}catch(const std::exception&e){std::cerr<<"ray_material error: "<<e.what()<<'\n';return 2;}}
