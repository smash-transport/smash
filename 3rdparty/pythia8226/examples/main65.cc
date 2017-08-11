// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program demonstrating the functionality of
// UserHooks which assign vertex information to MPI and shower particles.
// Author: Christian Bierlich.

#include "Pythia8/Pythia.h"
#include "Pythia8/UserHooks.h"

using namespace Pythia8;

//==========================================================================

// The UserHooks derived class
// This is only intended to serve as an example for how to use the vertex
// methods, not as a physical model for assigning production vertices.
class VertexUserHooks : public UserHooks{

public:

// The constructor takes parameters for the random dist of positions
VertexUserHooks(double MPIWidthIn, double emissionWidthIn)
  : widthMPI(MPIWidthIn), widthEmission(emissionWidthIn) {}

~VertexUserHooks() {}

// Must be set in order to enable the process
bool canSetProductionVertex() { return true; }

// Vertex for MPI production.
Vec4 vertexForMPI(Particle, double b) {
  bNow = b;
  pair<double, double> xy = rndmPtr->gauss2();
  Vec4 ret( bNow*widthMPI*xy.first, bNow*widthMPI*xy.second);
  return ret;
}

// Vertex for FSR production.
Vec4 vertexForFSR(Particle& p) {
  if (!p.hasVertex()) p.vProd(vertexForMPI(p,bNow));
  Vec4 v = p.vProd();
  Vec4 ret = getVertex(v);
  return ret;
}

// Vertex for ISR production.
Vec4 vertexForISR(Particle& p) {
  if (!p.hasVertex()) p.vProd(vertexForMPI(p,bNow));
  Vec4 v = p.vProd();
  Vec4 ret = getVertex(v);
  return ret;
}

// Smear vertex relative to mother.
Vec4 getVertex(Vec4 mother) {
  pair<double, double> xy = rndmPtr->gauss2();
  Vec4 ret( mother.px() + widthEmission * xy.first,
            mother.py() + widthEmission * xy.second);
  return ret;
}

private:
  double widthMPI, widthEmission, bNow;

};

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Create and give in UserHook to set vertex information.
  VertexUserHooks* userHooksPtr = new VertexUserHooks(0.2, 0.1);
  pythia.setUserHooksPtr( userHooksPtr);

  // Process selection. LHC initialization.
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("HadronLevel:all = off");
  pythia.init();

  // Begin event loop. Generate event, stop before hadronization.
  for (int iEvent = 0; iEvent < 1; ++iEvent) {
    if (!pythia.next()) continue;

    // Print vertex info.
    pythia.event.list();
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].vProd().px() == 0.0 ) cout << "   ";
      else cout << "-->" ;
      cout << pythia.event[i].index() << " " << pythia.event[i].status()
           << endl;
    }

  }

  // Finish.
  pythia.stat();
  return 0;
}
