#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace QuantLib;
using namespace std;


struct HestonPParams {
    Real S0;        // initial spot
    Real v0;        // initial variance
    Real mu;        // physical drift (P)

    Real kappaP;    // mean reversion speed (P)
    Real thetaP;    // long-run variance (P)
    Real xi;        // vol-of-vol (same in P/Q in your first model)
    Real rho;       // correlation (same in P/Q initially)
};

struct VRPParams {
    Real lambda;    // simplest: constant VRP for variance drift
};

struct HestonQParams {
    // State for pricing "today"
    Real S0;        // spot used for pricing (often S_t)
    Real v0;        // initial variance used for pricing (often v_t)

    // Market (Q drift uses r-q)
    Rate r;         // risk-free rate (or a curve later)
    Rate q;         // dividend yield (or carry)
    DayCounter dc;  // day count for time conversion

    // Q-model params
    Real kappaQ;
    Real thetaQ;
    Real xi;
    Real rho;
};

struct PPath{
    vector<Real> logS;
    vector<Real> v;
    vector<Real> returns;
};

struct CallGrid{
    vector <Date> maturities;
    vector <int> strikes;
    vector<vector<Real>> C;
    Real S0;
    Real v0;
};

void print_1d(vector<Real> arr){
    int n = arr.size();
    for (int i = 0; i < n; i++){
        cout << arr[i] <<",";
    }
}

ostream& operator<<(ostream& os, const PPath& p) {
    os << "\n=========== PPath ===========\n";
    os << "Path length      : " << p.logS.size() << "\n";
    os << "Variance length  : " << p.v.size() << "\n";
    os << "Returns length   : " << p.returns.size() << "\n\n";

    auto print_sample = [&](const std::vector<Real>& vec,
                            const std::string& name) {
        os << name << " (first 5): ";
        Size n = std::min<Size>(5, vec.size());
        for (Size i = 0; i < n; ++i)
            os << std::fixed << std::setprecision(6) << vec[i] << " ";
        if (vec.size() > 5) os << "...";
        os << "\n";
    };

    print_sample(p.logS, "logS");
    print_sample(p.v, "v");
    print_sample(p.returns, "returns");

    os << "================================\n";
    return os;
}


void print_Callgrid(const CallGrid& grid){

    cout <<"Starting S:"<<grid.S0<< endl;
    cout <<"Starting Vol:"<<grid.v0<< endl;
    
    cout << "printing inside the call grid"<<endl;
    cout<<"Call Price grid (rows = maturity, cols = strikes)"<<endl;
    cout << setw(20) << "T\\K";
    for (Real K : grid.strikes) cout << setw(12) << K;
    cout << "\n";

    for (Size iT = 0; iT < grid.maturities.size(); ++iT) {
        cout << setw(12) << grid.maturities[iT];
        for (Size iK = 0; iK < grid.strikes.size(); ++iK) {
            cout << setw(12) << fixed << setprecision(4) << grid.C[iT][iK];
        }
        cout << "\n";
    }

}

void print_Surfaces(const vector<CallGrid>& surfaces) {

    std::cout << "\n========== Heston Call Surfaces ==========\n";
    std::cout << "Total surfaces: " << surfaces.size() << "\n\n";

    for (Size i = 0; i < 5; ++i) {
        std::cout << "------------------------------------------\n";
        std::cout << "Surface ID: " << i << "\n";
        std::cout << "------------------------------------------\n";
        print_Callgrid(surfaces[i]);
        std::cout << "\n";
    }
}

ostream& operator<<(ostream& os, const HestonPParams& P) {
    os << "=== Heston P Parameters ===\n";
    os << "S0      : " << P.S0 << "\n";
    os << "v0      : " << P.v0 << "\n";
    os << "mu      : " << P.mu << "\n";
    os << "kappaP  : " << P.kappaP << "\n";
    os << "thetaP  : " << P.thetaP << "\n";
    os << "xi      : " << P.xi << "\n";
    os << "rho     : " << P.rho << "\n";
    os << "===========================\n";
    return os;
}

ostream& operator<<(ostream& os, const HestonQParams& Q) {
    os << "=== Heston Q Parameters ===\n";
    os << "S0      : " << Q.S0 << "\n";
    os << "v0      : " << Q.v0 << "\n";
    os << "r       : " << Q.r << "\n";
    os << "q       : " << Q.q << "\n";
    os << "DayCount: " << Q.dc.name() << "\n";
    os << "kappaQ  : " << Q.kappaQ << "\n";
    os << "thetaQ  : " << Q.thetaQ << "\n";
    os << "xi      : " << Q.xi << "\n";
    os << "rho     : " << Q.rho << "\n";
    os << "===========================\n";
    return os;
}


HestonQParams toQ(const HestonPParams& P, const VRPParams& V, Rate r, Rate q, const DayCounter& dc,
                const Real* S_today = nullptr, const Real* v_today = nullptr){
                    // Basic validation (keep it light but protective)
                    if (P.kappaP <= 0.0) throw std::invalid_argument("kappaP must be > 0");
                    if (P.thetaP <= 0.0) throw std::invalid_argument("thetaP must be > 0");
                    if (P.xi < 0.0)      throw std::invalid_argument("xi must be >= 0");
                    if (fabs(P.rho) > 1.0) throw std::invalid_argument("rho must be in [-1,1]");

                    Real kappaQ = P.kappaP + V.lambda;
                    if (kappaQ <= 0.0) throw invalid_argument("kappaQ must be > 0; Lambda is too negative?");
                    Real thetaQ = (P.kappaP * P.thetaP) / kappaQ;

                    HestonQParams Q;
                    Q.r = r;
                    Q.q = q;
                    Q.dc = dc;

                    Q.kappaQ = kappaQ;
                    Q.thetaQ = thetaQ;
                    Q.xi = P.xi;
                    Q.rho = P.rho;

                    Q.S0 = (S_today ? *S_today : P.S0);
                    Q.v0 = (v_today ? *v_today : P.v0);

                    if(Q.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
                    if(Q.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");
                    return Q;

}

void callgrid(HestonQParams& Q,CallGrid& grid){
    if (Q.kappaQ <= 0.0) throw invalid_argument("kappaP must be > 0");
    if (Q.thetaQ <= 0.0) throw invalid_argument("thetaP must be > 0");
    if (Q.xi < 0.0)      throw invalid_argument("xi must be >= 0");
    if (fabs(Q.rho) > 1.0) throw invalid_argument("rho must be in [-1,1]");
    if(grid.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
    if(grid.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");

    Calendar cal = TARGET();
    Date today(25,February,2026);
    Settings::instance().evaluationDate() = today;
    
    Handle<Quote> spot(ext::make_shared<SimpleQuote>(grid.S0));
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today, Q.r, Q.dc));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, Q.q, Q.dc));

    auto process = ext::make_shared<HestonProcess>(rTS,qTS,spot,grid.v0,Q.kappaQ,Q.thetaQ,Q.xi,Q.rho,HestonProcess::QuadraticExponential);
    auto model = ext::make_shared<HestonModel>(process);
    auto engine =  ext::make_shared<AnalyticHestonEngine>(model);

    //vector <vector <Real>> C(grid.maturities.size(),vector<Real> (grid.strikes.size(),0));
    grid.C.assign(grid.maturities.size(),vector<Real> (grid.strikes.size(),0));
    for (Size iT = 0; iT < grid.maturities.size();iT++){
        auto exercise = ext::make_shared<EuropeanExercise>(grid.maturities[iT]);
        for (Size iK = 0; iK < grid.strikes.size();iK ++ ){
            auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Call,grid.strikes[iK]);
            VanillaOption opt(payoff,exercise);
            opt.setPricingEngine(engine);
            grid.C[iT][iK] = opt.NPV();
        }
    }
    // cout<< "Printing inside the call surface function" <<endl;
    // cout<<"Call Price grid (rows = maturity, cols = strikes)"<<endl;
    // cout << setw(20) << "T\\K";
    // for (Real K : grid.strikes) cout << setw(12) << K;
    // cout << "\n";

    // for (Size iT = 0; iT < grid.maturities.size(); ++iT) {
    //     cout << setw(12) << grid.maturities[iT];
    //     for (Size iK = 0; iK < grid.strikes.size(); ++iK) {
    //         cout << setw(12) << fixed << setprecision(4) << C[iT][iK];
    //     }
    //     cout << "\n";
    // }

}


PPath logReturns(const HestonPParams& P, Size steps, Size seed){
    PPath ppath;
    if (P.kappaP <= 0.0) throw invalid_argument("kappaP must be > 0");
    if (P.thetaP <= 0.0) throw invalid_argument("thetaP must be > 0");
    if (P.xi < 0.0)      throw invalid_argument("xi must be >= 0");
    if (fabs(P.rho) > 1.0) throw invalid_argument("rho must be in [-1,1]");
    if(P.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
    if(P.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");

    Calendar cal = TARGET();
    Date today(26,February,2026);
    Settings::instance().evaluationDate() = today;
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today,P.mu,Actual365Fixed()));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, 0.0, Actual365Fixed()));
    Handle<Quote> s0(ext::make_shared<SimpleQuote>(P.S0));
    auto process = ext::make_shared<HestonProcess>(rTS,qTS,s0,P.v0,P.kappaP,P.thetaP,P.xi,P.rho,
                    HestonProcess::QuadraticExponential);

    Time T = 1.0;
    Time dt = T/steps;
    Size factors = 2;

    TimeGrid grid(T,steps);
    typedef PseudoRandom::rsg_type rsg_type;
    rsg_type rsg = PseudoRandom::make_sequence_generator(factors * grid.size()-factors, seed);
    MultiPathGenerator<rsg_type> gen(process,grid,rsg,false);
    const MultiPath path = gen.next().value;
    
    ppath.logS.resize(grid.size());
    ppath.v.resize(grid.size());
    ppath.returns.resize(grid.size() - 1);
    
    for(Size i = 0; i < grid.size();i++){
        ppath.logS[i] = log(path[0][i]);
        ppath.v[i] = path[1][i];
    }

    for(Size i = 1; i < grid.size();i++){
        ppath.returns[i-1] = ppath.logS[i] - ppath.logS[i-1]; 
    }
    return ppath;
}

int main(){
    HestonPParams P;
    HestonQParams Q;
    VRPParams VRP;
    //CallGrid grid;
    
    Date today(26, February, 2026);
    Calendar cal = TARGET();

    P.S0 = 100.0;
    P.v0 = 0.04;
    P.mu    = 0.05;   // physical drift you want
    P.kappaP = 1.5;
    P.thetaP = 0.04;
    P.xi    = 0.3;    // "sigma" of variance process in QL
    P.rho   = -0.7;
    Size steps = 252;
    Size seed = 42;

    PPath ppath = logReturns(P,steps,seed);
    
    VRP.lambda = 0.5;
    Rate r  = 0.03;   // risk-free (Q)
    Rate q  = 0.00;   // dividend yield (Q)
    const DayCounter dc = Actual365Fixed();
    // toQ(const HestonPParams& P, const VRPParams& V, Rate r, Rate q, const DayCounter& dc,
    //     const Real* S_today = nullptr, const Real* v_today = nullptr)
      
    Q = toQ(P,VRP,r,q,dc,&P.S0,&P.v0);
    cout << P <<endl;
    cout << Q <<endl;    
    cout << ppath <<endl;

    vector<CallGrid> surfaces;
    int SurfaceFrequency = 10;
    for(int i = 0; i < ppath.logS.size();i = i+SurfaceFrequency){
        CallGrid grid;
        grid.maturities = {
            cal.advance(today,1,Months),
            cal.advance(today,3,Months),
            cal.advance(today,6,Months),
            cal.advance(today,12,Months)
        };
    
        grid.strikes = {80,90,100,110,120};
        
        Real S_today = exp(ppath.logS[i]);
        Real v_today = ppath.v[i];
        
        grid.S0 = exp(ppath.logS[i]);
        grid.v0 = ppath.v[i];

        callgrid(Q,grid);
        surfaces.push_back(grid);
    }
    cout<<"Surfaces Shape:"<<surfaces.size()<<endl;
    print_Surfaces(surfaces);
        
}


// int main(){
//     HestonPParams P;
//     HestonQParams Q;
//     //CallGrid grid;
    
//     Date today(26, February, 2026);
//     Calendar cal = TARGET();

//     P.S0 = 100.0;
//     P.v0 = 0.04;
//     P.mu    = 0.05;   // physical drift you want
//     P.kappaP = 1.5;
//     P.thetaP = 0.04;
//     P.xi    = 0.3;    // "sigma" of variance process in QL
//     P.rho   = -0.7;
//     Size steps = 252;
//     Size seed = 42;
//     PPath ppath = logReturns(P,steps,seed);

//     // cout<< "Checking the paths inside the structure." <<endl;
//     // cout<<"log S:"<< endl;
//     // print_1d(ppath.logS);
//     // cout<<"v:"<< endl;
//     // print_1d(ppath.v);
//     // cout<<"log Returns:"<<endl;
//     // print_1d(ppath.returns);

//     Q.S0 = 100.0;
//     Q.r  = 0.03;   // risk-free (Q)
//     Q.q  = 0.00;   // dividend yield (Q)
//     Q.dc = Actual365Fixed();
//     Q.v0     = 0.04;
//     Q.kappaQ = 1.5;
//     Q.thetaQ = 0.04;
//     Q.xi     = 0.3;
//     Q.rho    = -0.7;
    
//     vector<CallGrid> surfaces;

//     // Real S0 = 100.0;
//     // Real v0 = 0.04; 
//     // grid.S0 = S0;
//     // grid.v0 = v0;

//     int SurfaceFrequency = 10;

//     for(int i = 0; i < ppath.logS.size();i = i+SurfaceFrequency){
//         CallGrid grid;
//         grid.maturities = {
//             cal.advance(today,1,Months),
//             cal.advance(today,3,Months),
//             cal.advance(today,6,Months),
//             cal.advance(today,12,Months)
//         };
    
//         grid.strikes = {80,90,100,110,120};
//         Real S_today = ppath.logS[i];
//         Real v_today = ppath.v[i];
//         callgrid(Q,grid,S_today,v_today);
//         surfaces.push_back(grid);
//     }
//     cout<<"Surfaces Shape:"<<surfaces.size()<<endl;


//     // grid.maturities = {
//     //     cal.advance(today,1,Months),
//     //     cal.advance(today,3,Months),
//     //     cal.advance(today,6,Months),
//     //     cal.advance(today,12,Months)
//     // };

//     // grid.strikes = {80,90,100,110,120};
//     // callgrid(Q,grid,S0,v0);
//     // print_Callgrid(grid);

// }

