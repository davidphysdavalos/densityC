//Plantilla lista pa chambear con itpp mas otras cosas
#include <iostream>
#include <cpp/dev_random.cpp>
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
#include <itpp/stat/histogram.h>
#include <cpp/itpp_ext_math.cpp>
#include <cpp/spinchain.cpp>
#include <itpp/stat/misc_stat.h>
#include <fstream>
#include <cpp/RMT.cpp>
#include <itpp/base/vec.h>

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
using namespace RMT;


TCLAP::CmdLine cmd("Programilla para calcular concurrencias", ' ', "0.9");
//TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"normalito", "string",cmd); // Para llamar strings
TCLAP::ValueArg<double> belltheta("","theta", "Tipo de estado separable, intermedio o bell segun sea theta" ,false, 0.6,"double",cmd);
TCLAP::ValueArg<double> purity("p","purity", "Purity",false, 0.6,"double",cmd);
TCLAP::ValueArg<string> model("","model", "Choose model" ,false,"total","string",cmd);
TCLAP::ValueArg<unsigned int> seed("","seed", "Random seed [0 for urandom]",false, 0,"unsigned int",cmd);
//TCLAP::ValueArg<int> members("m","members", "Ensemble members",false, 1,"int",cmd);
TCLAP::ValueArg<int> envdim("m","envdim", "Environment dimension",false, 4,"int",cmd);
TCLAP::ValueArg<int> spectatordim("","espdim", "Dimension of the spectator system",false, 4,"int",cmd);
TCLAP::ValueArg<int> realizations("r","realit", "Realizations",false, 100,"int",cmd);
TCLAP::ValueArg<double> epsilon("e","epsilon", "delta del tiempo",false, 0.001,"double",cmd);
TCLAP::ValueArg<double> coupling("c","coupling", "Coupling for spectator and so on",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> delta("d","delta", "delta del tiempo",false, 0.05,"double",cmd);
TCLAP::ValueArg<double> tiempo("t","tiempo", "tiempo final",false, 2.0,"double",cmd);
TCLAP::ValueArg<double> tiempocero("","t0", "tiempo inicial",false, 0.0,"double",cmd);
//TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd); //llamando enteros


int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

// {{{ Set seed for random
unsigned int semilla=seed.getValue();
if (semilla == 0){
  Random semilla_uran; semilla=semilla_uran.strong();
} 
RNG_reset(semilla);
// }}}


// state selection
cvec init=to_cvec(BellState(belltheta.getValue()));


cmat H;
cvec list [realizations.getValue()];
int totrace=(int) 3*envdim.getValue();

// Total Hamiltonian
if(model.getValue()=="total"){
cmat U;
int i;
cvec state;
cvec initextend=TensorProduct(init,RandomState(envdim.getValue()));
cmat partialstate;
cvec eigenval;
double pur;
int j;
for(int num=0;num<realizations.getValue();num++){
j=0;
i=1;
H=RandomGUE(4*envdim.getValue());
while(i==1){
U=exponentiate_nonsym(-complex <double>(0,1)*delta.getValue()*j*H);
j++;
state=U*initextend;
partialstate=partial_trace_qubits(state,totrace);
pur=(double)Purity(partialstate);
if(pur<=purity.getValue()+epsilon.getValue() and pur>=purity.getValue()-epsilon.getValue()){
	i=0;
}
}
eigenval=eig(partialstate);
cout<<Chop(real(eigenval(0)))<<' '<<Chop(real(eigenval(1)))<<' '<<Chop(real(eigenval(2)))<<' '<<Chop(real(eigenval(3)))<<endl;
}
}
// Spectator Hamiltonian

if(model.getValue()=="spectator"){
cmat U;
int i;
cvec state;
cvec initextend=TensorProduct(init,RandomState(envdim.getValue()));
cmat partialstate;
cvec eigenval;
double pur;
int j;
for(int num=0;num<realizations.getValue();num++){
j=0;
i=1;
H=TensorProduct(eye_c(4),RandomGUE(envdim.getValue()))+coupling.getValue()*TensorProduct(eye_c(2),RandomGUE(2*envdim.getValue()));
while(i==1){
U=exponentiate_nonsym(-complex <double>(0,1)*delta.getValue()*j*H);
j++;
state=U*initextend;
partialstate=partial_trace_qubits(state,totrace);
pur=(double)Purity(partialstate);
if(pur<=purity.getValue()+epsilon.getValue() and pur>=purity.getValue()-epsilon.getValue()){
	i=0;
}
}
//list[num]=eig(partialstate);
eigenval=eig(partialstate);
cout<<Chop(real(eigenval(0)))<<' '<<Chop(real(eigenval(1)))<<' '<<Chop(real(eigenval(2)))<<' '<<Chop(real(eigenval(3)))<<endl;
}
//}
//cvec eigenval;
//for(int num=0;num<realizations.getValue();num++){
//eigenval=list[num];
//cout<<Chop(real(eigenval(0)))<<' '<<Chop(real(eigenval(1)))<<' '<<Chop(real(eigenval(2)))<<' '<<Chop(real(eigenval(3)))<<endl;
//}
}

//Spectator 2 : Generalization to 3 qubits. Use GHZ and W to see what happens.

//Tuneable Coupling

if(model.getValue()=="tuneable"){
cmat U;
int i;
cvec state;
cvec initextend=TensorProduct(init,RandomState(envdim.getValue()));
cmat partialstate;
cvec eigenval;
double pur;
int j;
for(int num=0;num<realizations.getValue();num++){
j=0;
i=1;
H=TensorProduct(eye_c(4),RandomGUE(envdim.getValue()))+coupling.getValue()*RandomGUE(4*envdim.getValue());
while(i==1){
U=exponentiate_nonsym(-complex <double>(0,1)*delta.getValue()*j*H);
j++;
state=U*initextend;
partialstate=partial_trace_qubits(state,totrace);
pur=(double)Purity(partialstate);
if(pur<=purity.getValue()+epsilon.getValue() and pur>=purity.getValue()-epsilon.getValue()){
	i=0;
}
}
eigenval=eig(partialstate);
cout<<Chop(real(eigenval(0)))<<' '<<Chop(real(eigenval(1)))<<' '<<Chop(real(eigenval(2)))<<' '<<Chop(real(eigenval(3)))<<endl;
}
}


}
