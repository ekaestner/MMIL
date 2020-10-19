#include <iostream.h>
#include "Vertex.h"
#include "split.h"
#include "coplanar.h"
#include "outputDebug.h"


#include "TestCase.h"
#include "TestSuite.h"
#include "TestCaller.h"
#include "TestCaseRegister.h"

// 00=none 01=bold 04=underscore 05=blink 07=reverse 08=concealed
// 30=black 31=red 32=green 33=yellow 34=blue 35=magenta 36=cyan 37=white


inline void testIt( const char* const test_desc, bool test_case) 
{
  std::cout << test_desc;
  if (test_case) 
  {
    std::cout << "\t ["<<char(27)<<"[00;32m"<<"ok"<<char(27)<<"[00;00m"<<"]\n"; 
  }
  else 
  { 
     std::cout << "\t ["<<char(27)<<"[00;31m"<<"FAILED"<<char(27)<<"[00;00m"<<"]\n";
  }
}

class MathTestSuite : public TestCase
{
public:
   MathTestSuite( std::string name ) : TestCase( name )
   {
   }
   virtual ~MathTestSuite() {}
   virtual void setUp()
   {
   }

   static Test* newSuite()
   {
      TestSuite* testSuite = new TestSuite( "MathTestSuite" );

      testSuite->addTest( new TestCaller<MathTestSuite>( "detectAndFixCollinearVertices", detectAndFixCollinearVertices ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "detectAndFixIdenticalVerts", detectAndFixIdenticalVerts ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "split", split ));
      
      testSuite->addTest( new TestCaller<MathTestSuite>( "test10", test10 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test11", test11 ));
      
      testSuite->addTest( new TestCaller<MathTestSuite>( "test1", test1 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test1a", test1a ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test1b", test1b ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test1c", test1c ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test1d", test1d ));
      
      testSuite->addTest( new TestCaller<MathTestSuite>( "test2", test2 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test3", test3 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test4", test4 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test5", test5 ));
//      testSuite->addTest( new TestCaller<MathTestSuite>( "test6", test6 ));
      
      testSuite->addTest( new TestCaller<MathTestSuite>( "test7", test7 ));
      testSuite->addTest( new TestCaller<MathTestSuite>( "test9", test9 ));

      return testSuite;
   }


void detectAndFixCollinearVertices()
{
   Polygon<float> poly;
   poly.push_back( Vertex<float>( 0.0f, 0.0f, 0.0f ) );
   poly.push_back( Vertex<float>( 0.0f, 1.0f, 0.0f ) );
   poly.push_back( Vertex<float>( 0.0f, 2.0f, 0.0f ) );
   poly.push_back( Vertex<float>( 1.0f, 1.0f, 0.0f ) );
   poly.push_back( Vertex<float>( 0.0f, -2.0f, 0.0f ) );
   poly.push_back( Vertex<float>( 0.0f, -1.0f, 0.0f ) );

   testIt( "Polygon<type>::detectCollinear()", poly.detectCollinear() == true );
   assert( poly.detectCollinear() == true && "couldn't detect collinear verts" );

   poly.fixCollinear();
   
   testIt( "Polygon<type>::fixCollinear()", poly.detectCollinear() == false );
   assert( poly.detectCollinear() == false && "couldn't fix collinear verts" );
}

void detectAndFixIdenticalVerts()
{
{
   Polygon<float> polygon;
   polygon.push_back( Vertex<float>(91.7431f, 0.0f, 4.11116f) );
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );

   testIt( "Polygon<type>::detectIdenticalVerts() 1", polygon.detectIdenticalVerts() );
   assert( polygon.detectIdenticalVerts() == true && "should have detected the two identical verts" );
}
{
   Polygon<float> polygon;
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );
   polygon.push_back( Vertex<float>(91.7431f, 0.0f, 4.11116f) );
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );

   testIt( "Polygon<type>::detectIdenticalVerts() 2", polygon.detectIdenticalVerts() );
   assert( polygon.detectIdenticalVerts() == true && "should have detected the two identical verts" );
}
{
   Polygon<float> polygon;
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );
   polygon.push_back( Vertex<float>(91.7431f, 0.0f, 4.11116f) );

   testIt( "Polygon<type>::detectIdenticalVerts() 3", polygon.detectIdenticalVerts() );
   assert( polygon.detectIdenticalVerts() == true && "should have detected the two identical verts" );
}
}

void split()
{
   Polygon<float> polygon, positive, negative;
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 4.11116f) );
   polygon.push_back( Vertex<float>(91.7431f, 0.0f, 4.11116f) );
   polygon.push_back( Vertex<float>(91.7431f, 0.0f, 95.4704f) );
   polygon.push_back( Vertex<float>(91.7431f, 100.0f, 95.4704f) );

   Plane<float> plane;
   plane.set( Vec3<float>( 0.0f, -7.6294e-008, -1.0f ), -95.4704f );

   planeSplitsPolygon( plane, polygon, positive, negative );

   assert( negative.detectIdenticalVerts() == false && "split shouldn't cause this to happen" );

   int x;
   for (x = 0; x < positive.size(); ++x)
      std::cout<<"[split] pos: "<<(positive[x]).coord()[0]<<" "<<(positive[x]).coord()[1]<<" "<<(positive[x]).coord()[2]<<"\n"<<flush;
   for (x = 0; x < negative.size(); ++x)
      std::cout<<"[split] neg: "<<(negative[x]).coord()[0]<<" "<<(negative[x]).coord()[1]<<" "<<(negative[x]).coord()[2]<<"\n"<<flush;
}

void test1()
{
   Seg3<float> seg;
   Plane<float> plane;
   seg.set( Vec3<float>( 0.1f, 0.1f, -5.0f ), Vec3<float>( 0.1f, 0.1f, 5.0f ) );
   plane.set( Vec3<float>( 0,0,1 ), Vec3<float>( 0,0,0 ) );
   
   Vec3<float> intersectedPoint;
   bool result = false;
   intersectSegPlane( seg, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1] no intersection\n"<<flush;

   assert( result );
}

void test1a()
{
   bool result = false;

// easy test.
{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 1,1,0 ), 0.1 );
   seg.set( Vec3<float>( 1.0f, 100.0f, 0.0f ), Vec3<float>( 1.0f, -100.0f, 0.0f ) );

   Vec3<float> intersectedPoint;
   intersectSegPlane( seg, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1a] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1a] no intersection\n"<<flush;

   if (result == false)
   {
   const Vec3<float>& v0 = seg.p0();
   const Vec3<float>& v1 = seg.p1();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<" v0:    "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<" v1:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<" plane: "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<" vI:    "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
   }

   assert( result && "should be true" );
}
// another test.
{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 1,1,0 ), 1.41f );
   seg.set( Vec3<float>( -1.0f, -2.0f, 0.0f ), Vec3<float>( 3.0f, 2.0f, 0.0f ) );

   Vec3<float> intersectedPoint;
   intersectSegPlane( seg, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1a] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1a] no intersection\n"<<flush;

   if (result == false)
   {
   const Vec3<float>& v0 = seg.p0();
   const Vec3<float>& v1 = seg.p1();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<" v0:    "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<" v1:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<" plane: "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<" vI:    "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
   }

   assert( result && "should be true" );
   assert( intersectedPoint.isEqual( Vec3<float>( 1.49702, 0.49702, 0), 0.001f ) );
}
}

void test1c()
{
   bool result = false;


{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 1,1,0 ), 0.1 );
   seg.set( Vec3<float>( 1.0f, -1.0f, 0.0f ), Vec3<float>( 1.0f, 1.0f, 0.0f ) );

   Vec3<float> intersectedPoint;
   intersectSegPlane( seg, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1c] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1c] no intersection\n"<<flush;

if (result == false)
{
   const Vec3<float>& v0 = seg.p0();
   const Vec3<float>& v1 = seg.p1();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<"v0:    "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<"v1:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<"plane: "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<"vI:    "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
}


   assert( result && "should be true" );
}
}

#include "convert.h"
void test1b()
{
   bool result = false;

{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 1,1,0 ), 0.1 );
   seg.set( Vec3<float>( 1.0f, -1.0f, 0.0f ), Vec3<float>( 1.0f, 1.0f, 0.0f ) );

Ray3<float> ray;
seg2ray( seg, ray );

   Vec3<float> intersectedPoint;
   intersectRayPlane( ray, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1b] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1b] no intersection\n"<<flush;

if (result == false)
{
   const Vec3<float>& v0 = ray.origin();
   const Vec3<float>& v1 = ray.direction();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<"origin: "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<"dir:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<"plane:  "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<"vI:     "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
}

   assert( result && "should be true" );
}

{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 1,1,0 ), 0.1 );
   seg.set( Vec3<float>( 1.0f, 1.0f, 0.0f ), Vec3<float>( 1.0f, -1.0f, 0.0f ) );

Ray3<float> ray;
seg2ray( seg, ray );

   Vec3<float> intersectedPoint;
   intersectRayPlane( ray, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1b] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1b] no intersection\n"<<flush;

if (result == false)
{
   const Vec3<float>& v0 = ray.origin();
   const Vec3<float>& v1 = ray.direction();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<"origin: "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<"dir:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<"plane:  "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<"vI:     "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
}

   assert( result && "should be true" );
}

{
   Seg3<float> seg;
   Plane<float> plane( Vec3<float>( 0,1,0 ), 0.1 );
   seg.set( Vec3<float>( 1.0f, -1.0f, 0.0f ), Vec3<float>( 1.0f, 1.0f, 0.0f ) );

Ray3<float> ray;
seg2ray( seg, ray );

   Vec3<float> intersectedPoint;
   intersectRayPlane( ray, plane, intersectedPoint, result );

   if (result)
      std::cout<<"[test1b] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test1b] no intersection\n"<<flush;

if (result == false)
{
   const Vec3<float>& v0 = ray.origin();
   const Vec3<float>& v1 = ray.direction();
   const Vec3<float>& vI = intersectedPoint;
   std::cout<<"v0:    "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<"\n"<<flush;
   std::cout<<"v1:    "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<"\n"<<flush;
   std::cout<<"plane: "<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;
   std::cout<<"vI:    "<<vI[0]<<" "<<vI[1]<<" "<<vI[2]<<"\n"<<flush;
}


   assert( result && "should be true" );
}
}

// coplanar
void test1d()
{
   Plane<float> plane1( Vec3<float>( 0.0f, 1.0f, 0.0f ), Vec3<float>( 0.0f, 100.0f, 0.0f ) );
   Plane<float> plane2( Vec3<float>( 1.0f, 1.0f, 0.0f ), Vec3<float>( 0.0f, 100.0f, 0.0f ) );
   Plane<float> plane3( Vec3<float>( 0.0f, 1.0f, 0.0f ), Vec3<float>( 0.0f, 99.0f, 0.0f ) );
   Plane<float> plane4( Vec3<float>( 0.0f, -1.0f, 0.0f ), Vec3<float>( 0.0f, 100.0f, 60.4f ) );
   Plane<float> plane5( Vec3<float>( 0.0f, 1.0f, 0.0f ), Vec3<float>( 100.0f, 100.0f, 0.0f ) );
   
   assert( isCoplanar( plane1, plane1 ) && "should be coplanar" );
   assert( isCoplanar( plane1, plane5 ) && "should be coplanar" );
   
   assert( !isCoplanar( plane1, plane2 ) && "should not be coplanar" );
   assert( !isCoplanar( plane1, plane3 ) && "should not be coplanar" );
   
   assert( isCoplanarOpposite( plane1, plane4 ) && "should be coplanar but opposite" );
   assert( isCoplanarOpposite( plane4, plane5 ) && "should be coplanar but opposite" );
}

void test2()
{
   Tri3<float> tri;
   tri.p0( 0, 0, 0 );
   tri.p1( 0, 2, 0 );
   tri.p2( 2, 0, 0 );

   Seg3<float> seg;
   seg.set( Vec3<float>( 0.1f, -0.1f, -5.0f ), Vec3<float>( 0.1f, 0.1f, 5.0f ) );

   bool result = false;
   Vec3<float> intersectedPoint;
   intersectSegTri( seg, tri, intersectedPoint, result ); /// no work, test me

   if (result)
      std::cout<<"[test2] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test2] no intersection\n"<<flush;

   
   
   tri.p0( 0, 0, 0 );
   tri.p1( 0, -2, 0 );
   tri.p2( 2, 0, 0 );
   

   kev::intersectSegTri( seg, tri, intersectedPoint, result ); /// no work, test me

   if (result)
      std::cout<<"[test2] Intersected: "<<intersectedPoint[0]<<" "<<intersectedPoint[1]<<" "<<intersectedPoint[2]<<"\n"<<flush;
   else
      std::cout<<"[test2] no intersection\n"<<flush;

assert( result );
}

void  test3()
{
   Plane<float> plane;
   plane.set( Vec3<float>( 1.0f, 1.0f, 0.0f ), Vec3<float>( 0.0f, 0.1f, 0.0f ) );

   Polygon<float> polygon, positive, negative;
   polygon.push_back( Vertex<float>( -1,-1,0 ) );
   polygon.push_back( Vertex<float>( 1,-1,0 ) );
   polygon.push_back( Vertex<float>( 1,1,0 ) );
   polygon.push_back( Vertex<float>( -1,1,0 ) );

   planeSplitsPolygon( plane, polygon, positive, negative );

   int x;
   for (x = 0; x < positive.size(); ++x)
      std::cout<<"[test3] pos: "<<(positive[x]).coord()[0]<<" "<<(positive[x]).coord()[1]<<" "<<(positive[x]).coord()[2]<<"\n"<<flush;
   for (x = 0; x < negative.size(); ++x)
      std::cout<<"[test3] neg: "<<(negative[x]).coord()[0]<<" "<<(negative[x]).coord()[1]<<" "<<(negative[x]).coord()[2]<<"\n"<<flush;
   
}



void testSplit( int testnum, Polygon<float>& polygon, const Polygon<float>& cutter, int possize, int negsize )
{
   Plane<float> plane;
   cutter.getPlane( plane );
   
   Polygon<float> positive, negative;
   planeSplitsPolygon( plane, polygon, positive, negative );
   
   std::cout<<"[test"<<testnum<<"] pos: "; positive.outputDebug(); std::cout<<"\n"<<flush;
   std::cout<<"[test"<<testnum<<"] neg: "; negative.outputDebug(); std::cout<<"\n"<<flush;

   assert( possize == positive.size() && negsize == negative.size() && "TEST FAILED: sizes are unexpected");
}

// cut a polygon that has one edge against the other (shouldn't cut it at all)
void test4()
{
   Polygon<float> polygon, cutter;
   
// NOTE: polys are CCW winding

   polygon.push_back( Vertex<float>( 0.5, 0.0, 0.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 1.0, 0.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 0.0, 5.0 ) );

   cutter.push_back( Vertex<float>( 0.0, 0.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 0.0, 1.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 1.0, 0.0, 0.0 ) );
   
   testSplit( 4, polygon, cutter, 0, 3 );
}

// colinear verts?
void test6()
{
   Polygon<float> polygon, cutter;
   
   polygon.push_back( Vertex<float>( 0.5, 0.0, 0.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 0.5, 0.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 1.0, 0.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 0.0, 5.0 ) );

   cutter.push_back( Vertex<float>( 0.0, 0.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 0.0, 1.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 1.0, 0.0, 0.0 ) );
   
   testSplit( 6, polygon, cutter, 0, 3 );
}

// cut a polygon (should give negsize=4, possize=3)
void test5()
{
   Polygon<float> polygon, cutter;
   
   polygon.push_back( Vertex<float>( 0.5, 0.0, -1.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 1.0, -1.0 ) );
   polygon.push_back( Vertex<float>( 0.5, 0.0,  5.0 ) );

   cutter.push_back( Vertex<float>( 0.0, 0.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 0.0, 1.0, 0.0 ) );
   cutter.push_back( Vertex<float>( 1.0, 0.0, 0.0 ) );
   
   testSplit( 5, polygon, cutter, 4, 3 );
}


void test7()
{
   int size = 5;
   for (int x=0; x < size * 3; ++x)
   {
      int x_minus_one, x_current, x_plus_one, x_plus_two;
      getVectorNeighbors( size, x, x_minus_one, x_current, x_plus_one, x_plus_two );
      std::cout<<x_minus_one<<" "<<x_current<<" "<<x_plus_one<<" "<<x_plus_two<<"\n"<<flush;
   }
}

void test9()
{
   Vec3<float> vec1( 0, 36, 0 ), vec2( 0,0,0 );
   Seg3<float> seg;
   seg.set( vec1, vec2 );

   Plane<float> plane; 
   plane.set( Vec3<float>( 0, 1 ,0 ), Vec3<float>(0, 36, 0) );

   Vec3<float> intersectedPoint;
   bool result;
   intersectSegPlane( seg, plane, intersectedPoint, result );
   assert( result == true && "seg should have intersected plane" );

   seg.set( vec2, vec1 );
   intersectSegPlane( seg, plane, intersectedPoint, result );
   assert( result == true && "order shouldn't matter" );

   std::cout<<"[test9] insection seg/plane worked when seg end points were on the plane\n"<<flush; 
}

void test10()
{
   {
      Plane<float> plane( Vec3<float>( 0.0f, 0.0f, 1.0f ), Vec3<float>( 0.0f, 0.0f, 0.0f ) );

      float d;
      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, 11.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 11.0f );

      d = plane.getDistance( Vec3<float>( 0.0f, 10.0f, -20.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -20.0f );

      d = plane.getDistance( Vec3<float>( 203.0f, 22.0f, 1000.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 1000.0f );

      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, -11.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -11.0f );
   }
   {
      Plane<float> plane( Vec3<float>( 1.0f, 1.0f, 0.0f ), 0.1 );

      float d;
      // test on the plane, dist should be 0
      d = plane.getDistance( Vec3<float>( 0.070710678f, 0.070710678f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 0.0f );

      // test from the origin, dist should be -0.1 (neg side of plane)
      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -0.1f );

      d = plane.getDistance( Vec3<float>( -0.070710678f, -0.070710678f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -0.2f );
   }
   {
      Plane<float> plane( Vec3<float>( 1.0f, 1.0f, 0.0f ), Vec3<float>( 0.070710678f, 0.070710678f, 0.0f ) );

      float d;
      // test on the plane, dist should be 0
      d = plane.getDistance( Vec3<float>( 0.070710678f, 0.070710678f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 0.0f );

      // test from the origin, dist should be -0.1 (neg side of plane)
      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -0.1f );

      d = plane.getDistance( Vec3<float>( -0.070710678f, -0.070710678f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == -0.2f );
   }
   
   {
      Plane<float> plane( Vec3<float>( 0.0f, 0.0f, 1.0f ), Vec3<float>( 0.0f, 0.0f, -1000.0f ) );

      float d;
      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, 11.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 1011.0f );

      d = plane.getDistance( Vec3<float>( 0.0f, 10.0f, -20.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 980.0f );

      d = plane.getDistance( Vec3<float>( 203.0f, 22.0f, 0.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 1000.0f );

      d = plane.getDistance( Vec3<float>( 0.0f, 0.0f, -11.0f ) );
      std::cout<<d<<"\n"<<flush;
      assert( d == 989.0f );
   }
   
   
   {
      Plane<float> plane( Vec3<float>( 0.0f, 1.0f, 1.0f ), Vec3<float>( 0.0f, 0.0f, 0.0f ) );
      float d1, d2;
      d1 = plane.getDistance( Vec3<float>( 0.0f, 100.0f, 1000.0f ) );
      std::cout<<d1<<"\n"<<flush;
      
      d2 = plane.getDistance( Vec3<float>( 0.0f, -256.0f, 1000.0f ) );
      std::cout<<d2<<"\n"<<flush;
      
      assert( d1 > d2 );
   }
}

void test11()
{
   {
      Polygon<float> polygon;

      polygon.push_back( Vertex<float>( 0.0, 0.0, 100.0 ) );
      polygon.push_back( Vertex<float>( 1.0, 0.0, 100.0 ) );
      polygon.push_back( Vertex<float>( 1.0, 1.0, 100.0 ) );

      Vec3<float> norm;
      Plane<float> plane;
      polygon.getNormal( norm );
      polygon.getPlane( plane );

      std::cout<<norm[0]<<" "<<norm[1]<<" "<<norm[2]<<"\n"<<flush;
      std::cout<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;

      assert( norm == Vec3<float>( 0.0f, 0.0f, 1.0f ) );
      assert( plane.normal() == Vec3<float>( 0.0f, 0.0f, 1.0f ) );
      assert( plane.constant() == 100 );
   }
   
   {
      Polygon<float> polygon;

      polygon.push_back( Vertex<float>( -200.0, 0.0, 0.0 ) );
      polygon.push_back( Vertex<float>( -200.0, 1.0, 0.0 ) );
      polygon.push_back( Vertex<float>( -200.0, 0.0, 1.0 ) );

      Vec3<float> norm;
      Plane<float> plane;
      polygon.getNormal( norm );
      polygon.getPlane( plane );

      std::cout<<norm[0]<<" "<<norm[1]<<" "<<norm[2]<<"\n"<<flush;
      std::cout<<plane.normal()[0]<<" "<<plane.normal()[1]<<" "<<plane.normal()[2]<<" : "<<plane.constant()<<"\n"<<flush;

      assert( norm == Vec3<float>( 1.0f, 0.0f, 0.0f ) );
      assert( plane.normal() == Vec3<float>( 1.0f, 0.0f, 0.0f ) );
      assert( plane.constant() == -200 );
   }
}

};


TestCaseRegister<MathTestSuite> bokbok( "MathTestSuite" );



/*
void main()
{
   collinear();
   split();
   detect();

  test10();
  test11();


   test1a();
   test1b();
   test1c();
   test1d();

   test1();
   test2();
   test7();
   test3();
   test4();
   test5();
  // test6();
   test9();
   
}

*/
