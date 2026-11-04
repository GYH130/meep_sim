import com.comsol.model.*;
import com.comsol.model.physics.PhysicsFeature;
import com.comsol.model.util.ModelUtil;
import java.time.Instant;
import java.security.MessageDigest;
import java.util.*;

/** Independent explicit cavity-FLOOR-FILL pilot, NOT conformal wall/land coating.
 * Closed eleven-case H80 angular-preview/polarization candidate: six s
 * complements and five p-only preview-grid fills,3um floorfill,10.5/13um.
 * Actual s/p global vectors and mode axes require dedicated native admission.
 * No full hemisphere, complete-band, unpolarized angular-cone or accuracy claim.
 * Local source audit is not native compilation, source admission or qualification.
 * exact direct out-of-core48GiB configuration persisted in each mesh artifact.
 * Derived from frozen e6f202 CavityTiPilot; original source is unchanged.
 * A3 lower0..selected1/2/3um frustum is SiO2; upper cavity/cap air, original Ti surfaces.
 * A hash-bound manual admission is required before any physical runAll().
 * No Java filesystem, environment, property, network, or process access.
 * The controller embeds run(...) into a standalone CaseEntry.java before compiling.
 * All geometry numbers and mesh coordinates are SI metres (geom lengthUnit m).
 * References: official 6.3 RF Hexagonal Grating and RF Port documentation;
 * installed 6.3 Wave Optics hexagonal_grating dmodel confirms order property names.
 * This is a method-feasibility pilot, not angular/spectral emissivity qualification.
 */
public class FlatPrecheckEntry {
    public static long reviewedDof=-1L,reviewedTet=-1L;
    public static String reviewedMeshSha="";
    private static final double ABS_POWER_GATE=.005, REL_LOSS_GATE=.02, STRICT_ABS_DIAGNOSTIC=.002;
    private static final double UM=1e-6, PITCH=80*UM, APOTHEM=40*UM;
    private static final double AOPEN=69.5455842631299*UM/2, ABOTTOM=7.5*UM;
    private static final double TOL=1e-11, GEOMETRY_MEASURE_REL=1e-6;
    private static final String PROOF_SCHEMA="ANGLE_BAND_SOURCE_POWER_V2";
    private static final String INCIDENT_POWER_METHOD="FIXED_KZ_BLOCH_DEPHASED_FIVE_PLANE_TANGENTIAL_SP_COMPLEX_LS_V2";
    // Same-reference native hex probe measured alpha2=0 -> global phi=150deg.
    // The -150deg correction is additionally rechecked by the 30/60deg flat
    // controls and every production Variables-stage actual-k admission.
    public static final double PORT_AZIMUTH_OFFSET_DEG=-150;
    private static final double[] POWER_FIT_CAP_FRACTIONS={.1,.25,.5,.75,.9};
    private static final double POWER_FIT_CONDITION_MAX=10, POWER_FIT_RELATIVE_L2_MAX=.005, POWER_FIT_LOO_PIN_ERROR_MAX=.005;
    private static final double MODE_DIRECTION_LEAKAGE_MAX=1e-6;
    private static String root, caseId, started;
    private static int evalId=0;
    private static String lastPhase, lastExpression, lastNumericalTag, lastNumericalKind;
    private static boolean solverCompleted, solvedArtifactSaved;

    private static String json(Object o) {
        if(o==null)return "null";
        if(o instanceof String)return "\""+((String)o).replace("\\","\\\\").replace("\"","\\\"").replace("\n","\\n").replace("\r","\\r")+"\"";
        if(o instanceof Number) {double n=((Number)o).doubleValue();return Double.isFinite(n)?o.toString():"null";}
        if(o instanceof Boolean)return o.toString();
        if(o instanceof Map) {StringJoiner s=new StringJoiner(",","{","}");for(Object k:((Map<?,?>)o).keySet())s.add(json(k.toString())+":"+json(((Map<?,?>)o).get(k)));return s.toString();}
        if(o instanceof Collection) {StringJoiner s=new StringJoiner(",","[","]");for(Object x:(Collection<?>)o)s.add(json(x));return s.toString();}
        if(o instanceof int[]) {List<Integer> a=new ArrayList<>();for(int x:(int[])o)a.add(x);return json(a);}
        if(o instanceof double[]) {List<Double> a=new ArrayList<>();for(double x:(double[])o)a.add(x);return json(a);}
        throw new IllegalArgumentException("Unsupported JSON value "+o.getClass().getName());
    }
    private static Map<String,Object> obj(Object... kv) {Map<String,Object> m=new LinkedHashMap<>();for(int i=0;i<kv.length;i+=2)m.put((String)kv[i],kv[i+1]);return m;}
    private static void emit(String marker,Map<String,Object> m) {if(!m.containsKey("case_id"))m.put("case_id",caseId);System.out.println(marker+"="+json(m));System.out.flush();}
    private static void progress(String state) {lastPhase=state;emit("CODEX_PROGRESS_JSON",obj("case_id",caseId,"state",state,"started_utc",started,"updated_utc",Instant.now().toString()));}
    private static void evaluationCheckpoint(String phase,Map<String,Object> details) {
        lastPhase=phase;
        Map<String,Object> receipt=obj("case_id",caseId,"phase",phase,"solver_completed",solverCompleted,
            "solved_artifact_saved",solvedArtifactSaved,"qualification","UNVALIDATED_EVALUATION_CHECKPOINT",
            "numerical_evaluations_started",evalId,"last_numerical_tag",lastNumericalTag,
            "last_numerical_kind",lastNumericalKind,"last_expression",lastExpression,
            "updated_utc",Instant.now().toString());
        receipt.putAll(details);emit("CODEX_EVAL_CHECKPOINT_JSON",receipt);
    }
    private static int[] ints(Collection<Integer> a) {int[] r=new int[a.size()];int j=0;for(int n:a)r[j++]=n;return r;}
    private static int[] combine(int[]... lists) {TreeSet<Integer>s=new TreeSet<>();for(int[] a:lists)for(int x:a)s.add(x);return ints(s);}
    private static boolean close(double x,double y,double rel) {return Math.abs(x-y)<=Math.max(Math.abs(y)*rel,1e-24);}
    private static boolean allZ(double[][] v,int[] ids,double z) {for(int i:ids)if(Math.abs(v[2][i-1]-z)>TOL)return false;return ids.length>0;}
    private static boolean side(double[][] v,int[] ids,double ang,double zmin) {for(int i:ids)if(Math.abs(v[0][i-1]*Math.cos(ang)+v[1][i-1]*Math.sin(ang)-APOTHEM)>TOL||v[2][i-1]<zmin-TOL)return false;return ids.length>0;}
    private static double area(Model m,int[] faces) {if(faces.length==0)return 0;GeomMeasureFinal q=m.component("comp1").geom("geom1").measureFinal();q.selection().geom("geom1",2);q.selection().set(faces);return q.getArea();}
    private static void explicit(Model m,String tag,int dim,int[] ids) {m.component("comp1").selection().create(tag,"Explicit");m.component("comp1").selection(tag).geom("geom1",dim);m.component("comp1").selection(tag).set(ids);}
    private static void hex(GeomSequence g,String wp,String ext,double ap,double z,double distance,double scale) {
        g.create(wp,"WorkPlane");g.feature(wp).set("quickz",Double.toString(z));
        g.feature(wp).geom().create("hex","Polygon");
        double[] x=new double[6],y=new double[6];
        for(int j=0;j<6;j++){double t=Math.PI/6+j*Math.PI/3;x[j]=ap/Math.cos(Math.PI/6)*Math.cos(t);y[j]=ap/Math.cos(Math.PI/6)*Math.sin(t);}
        g.feature(wp).geom().feature("hex").set("x",x);g.feature(wp).geom().feature("hex").set("y",y);
        g.run(wp);g.create(ext,"Extrude");g.feature(ext).set("workplane",wp);
        g.feature(ext).set("distance",distance);g.feature(ext).set("scale",new double[][]{{scale,scale}});g.run(ext);
    }
    private static double integral(Model m,int[] faces,String expr) {
        if(faces.length==0)return 0;
        String t="ev"+(++evalId);m.result().numerical().create(t,"IntSurface");
        lastNumericalTag=t;lastNumericalKind="IntSurface";lastExpression=expr;
        m.result().numerical(t).selection().set(faces);m.result().numerical(t).set("data","dset1");
        m.result().numerical(t).set("expr",expr);m.result().numerical(t).set("intorderactive","on");m.result().numerical(t).set("intorder",6);
        double q=m.result().numerical(t).getReal()[0][0];m.result().numerical().remove(t);return q;
    }
    private static double avg(Model m,int[] top,double a,String expr) {return integral(m,top,expr)/a;}
    private static double sq(double re,double im){return re*re+im*im;}
    private static double global(Model m,String expr) {String t="ev"+(++evalId);lastNumericalTag=t;lastNumericalKind="EvalGlobal";lastExpression=expr;m.result().numerical().create(t,"EvalGlobal");m.result().numerical(t).set("data","dset1");m.result().numerical(t).set("expr",expr);double q=m.result().numerical(t).getReal()[0][0];m.result().numerical().remove(t);return q;}
    public static double[] requestedIncidentK(double lambdaUm,double thetaDeg,double phiDeg) {
        if(!Double.isFinite(lambdaUm)||lambdaUm<=0||!Double.isFinite(thetaDeg)||thetaDeg<0||thetaDeg>=90||!Double.isFinite(phiDeg))throw new IllegalArgumentException("Finite positive wavelength and nongrazing incidence required");
        double k=2*Math.PI/(lambdaUm*UM),t=Math.toRadians(thetaDeg),p=Math.toRadians(phiDeg);
        return new double[]{k*Math.sin(t)*Math.cos(p),k*Math.sin(t)*Math.sin(p),-k*Math.cos(t)};
    }
    public static double tangentialPlanePower(double area,double Z0,double cosTheta,double[] p,double[] s) {
        if(!Double.isFinite(area)||area<=0||!Double.isFinite(Z0)||Z0<=0||!Double.isFinite(cosTheta)||cosTheta<=0||cosTheta>1+1e-12)throw new IllegalArgumentException("Physical air-port admittance required");
        return area*(sq(p[0],p[1])/cosTheta+cosTheta*sq(s[0],s[1]))/(2*Z0);
    }
    public static void setProbeContext(String id){caseId=id;started=Instant.now().toString();evalId=0;}
    /** Pure arithmetic helper. Return [a.re,a.im,b.re,b.im,SSE,sum|E|^2,cond2(X)].
     * E(z)=a*exp(+i*k0*z)+b*exp(-i*k0*z); k0 is fixed, never fitted.
     * The normal-equation inverse is analytic Hermitian 2x2; condition is for
     * the design matrix X itself, not its squared Gram-matrix condition.
     */
    public static double[] fitComplexPlaneWaves(double[] z,double[][] e,double k0) {
        if(z==null||e==null||z.length!=e.length||z.length<2||!Double.isFinite(k0))throw new IllegalArgumentException("Finite matching plane samples required");
        int n=z.length;double cr=0,ci=0,pr=0,pi=0,qr=0,qi=0,energy=0;
        for(int j=0;j<n;j++){
            if(e[j]==null||e[j].length!=2||!Double.isFinite(z[j])||!Double.isFinite(e[j][0])||!Double.isFinite(e[j][1]))throw new IllegalArgumentException("Nonfinite complex plane sample");
            double c=Math.cos(k0*z[j]),s=Math.sin(k0*z[j]),er=e[j][0],ei=e[j][1];
            cr+=Math.cos(2*k0*z[j]);ci-=Math.sin(2*k0*z[j]);
            pr+=c*er+s*ei;pi+=c*ei-s*er;qr+=c*er-s*ei;qi+=c*ei+s*er;energy+=er*er+ei*ei;
        }
        double h=Math.hypot(cr,ci),low=n-h,high=n+h;
        if(low<=1e-12*n)throw new IllegalArgumentException("Singular or numerically rank-deficient two-wave design");
        double det=low*high;
        double ar=(n*pr-cr*qr+ci*qi)/det,ai=(n*pi-cr*qi-ci*qr)/det;
        double br=(n*qr-cr*pr-ci*pi)/det,bi=(n*qi-cr*pi+ci*pr)/det,sse=0;
        for(int j=0;j<n;j++){double c=Math.cos(k0*z[j]),s=Math.sin(k0*z[j]);double er=(ar+br)*c+(bi-ai)*s,ei=(ai+bi)*c+(ar-br)*s;sse+=sq(er-e[j][0],ei-e[j][1]);}
        return new double[]{ar,ai,br,bi,sse,energy,Math.sqrt(high/low)};
    }
    /** Solver-free arithmetic contract, also callable by the native flat probe.
     * Independent expected counts and genuine global s/p source are fixed for eleven IDs.
     * This does NOT replace the measured hex-port orientation/power calibration.
     */
    public static Map<String,Object> pureArithmeticSelfTest() {
        List<Map<String,Object>> rows=new ArrayList<>();boolean pass=true;
        String[] ids={"ML6_A3H100_t2_l10p5_th0_phi0_s","ML7_A3H100_t2_l10p5_th0_phi0_p","ML8_A3H100_t2_l10p5_th30_phi30_s","ML9_A3H100_t2_l10p5_th30_phi30_p"};
        String[] polarizations={"s","p","s","p"};
        double[] lambdas={10.5,10.5,10.5,10.5},angles={0,0,30,30},azimuths={0,0,30,30},heights={100,100,100,100};int[] counts={163,163,157,157},padded={199,199,190,190};
        for(int c=0;c<lambdas.length;c++){
            double k0=2*Math.PI/(lambdas[c]*UM),ct=Math.cos(Math.toRadians(angles[c]));double[] kin=requestedIncidentK(lambdas[c],angles[c],azimuths[c]),kpad=requestedIncidentK(lambdas[c]/1.1,angles[c],azimuths[c]);
            int n=expectedOrders(lambdas[c],kin[0],kin[1]).size(),npad=expectedOrders(lambdas[c]/1.1,kpad[0],kpad[1]).size();
            double[]z=new double[POWER_FIT_CAP_FRACTIONS.length];double[][]e=new double[z.length][2];
            for(int j=0;j<z.length;j++){z[j]=heights[c]*UM+POWER_FIT_CAP_FRACTIONS[j]*10.5*UM;double q=k0*ct*z[j],co=Math.cos(q),si=Math.sin(q);e[j][0]=1.3*co-.3*si;e[j][1]=.1*co+.7*si;}
            double[]f=fitComplexPlaneWaves(z,e,k0*ct);double error=Math.max(Math.max(Math.abs(f[0]-1),Math.abs(f[1]-.2)),Math.max(Math.abs(f[2]-.3),Math.abs(f[3]+.1)));
            double a=2*Math.sqrt(3)*APOTHEM*APOTHEM,Z0=376.730313668,pamp=Math.sqrt(2*Z0*ct/a),samp=Math.sqrt(2*Z0/(a*ct));
            double pp=tangentialPlanePower(a,Z0,ct,new double[]{pamp,0},new double[]{0,0}),ps=tangentialPlanePower(a,Z0,ct,new double[]{0,0},new double[]{samp,0});
            double[] targetE=globalE(polarizations[c],angles[c],azimuths[c],false);double transverseDot=kin[0]*targetE[0]+kin[1]*targetE[1]+kin[2]*targetE[2];
            double phiRad=Math.toRadians(azimuths[c]),thetaRad=Math.toRadians(angles[c]);double[] independentTarget="s".equals(polarizations[c])?new double[]{-Math.sin(phiRad),Math.cos(phiRad),0}:new double[]{Math.cos(thetaRad)*Math.cos(phiRad),Math.cos(thetaRad)*Math.sin(phiRad),Math.sin(thetaRad)};double sourceVectorError=0;for(int j=0;j<3;j++)sourceVectorError=Math.max(sourceVectorError,Math.abs(targetE[j]-independentTarget[j]));
            boolean ok=authorizedPolarization(ids[c]).equals(polarizations[c])&&authorizedFilmUm(ids[c])==new int[]{2,2,2,2}[c]&&authorizedHeightUm(ids[c])==heights[c]&&authorizedLambda(ids[c])==lambdas[c]&&authorizedTheta(ids[c])==angles[c]&&authorizedPhi(ids[c])==azimuths[c]&&authorizedPaddingChannels(ids[c])==2*padded[c]&&n==counts[c]&&npad==padded[c]&&error<1e-12&&f[4]<1e-20&&Math.abs(pp-1)<1e-12&&Math.abs(ps-1)<1e-12&&Math.abs(transverseDot)<1e-9&&sourceVectorError<1e-14;
            rows.add(obj("case_id",ids[c],"polarization",polarizations[c],"target_incoming_E",targetE,"height_um",heights[c],"source_plane_z_m",z,"lambda_um",lambdas[c],"theta_deg",angles[c],"phi_global_deg",azimuths[c],"physical_spatial_orders",n,"padded_spatial_orders",npad,"fit_max_coefficient_error",error,"fit_SSE",f[4],"p_power_W",pp,"s_power_W",ps,"incoming_k_dot_E",transverseDot,"independent_source_vector_max_abs_error",sourceVectorError,"passed",ok));pass&=ok;
        }
        boolean singularRejected=false;try{fitComplexPlaneWaves(new double[]{0,0},new double[][]{{1,0},{1,0}},1);}catch(IllegalArgumentException e){singularRejected=true;}
        boolean closedIdsReject=false;try{authorizedFilmUm("W1_A3_t2_p");}catch(IllegalArgumentException e){closedIdsReject=true;}
        pass&=singularRejected&&closedIdsReject;
        Map<String,Object> proof=obj("proof_schema",PROOF_SCHEMA,"state",pass?"PASS_PURE_ARITHMETIC_SELF_TEST":"FAIL_PURE_ARITHMETIC_SELF_TEST","passed",pass,"rows",rows,"singular_design_rejected",singularRejected,"unapproved_case_rejected",closedIdsReject,"solver_executed",false,"does_not_establish_native_port_orientation_or_cavity_acceptance",true);
        if(!pass)throw new IllegalStateException("Angle-band pure arithmetic self-test failed: "+json(proof));return proof;
    }
    private static double[] fittedValue(double[] fit,double k0,double z){double c=Math.cos(k0*z),s=Math.sin(k0*z);return new double[]{(fit[0]+fit[2])*c+(fit[3]-fit[1])*s,(fit[1]+fit[3])*c+(fit[0]-fit[2])*s};}
    private static double[] complexCutIntegral(Model m,String data,String expr) {
        String t="ev"+(++evalId);lastNumericalTag=t;lastNumericalKind="IntSurfaceCutPlane";lastExpression=expr;
        m.result().numerical().create(t,"IntSurface");
        try{NumericalFeature q=m.result().numerical(t);q.set("data",data);q.selection().all();q.set("expr",expr);q.set("intorderactive","on");q.set("intorder",8);double[][][]v=q.getComplex();return new double[]{v[0][0][0],v.length>1&&v[1]!=null?v[1][0][0]:0};}
        finally{m.result().numerical().remove(t);}
    }
    public static Map<String,Object> internalElectricPowerFit(Model m,double area,double zland,double cap,double lambda,double[] incidentK) {
        double k0=2*Math.PI/(lambda*UM),z0=global(m,"Z0_const");int n=POWER_FIT_CAP_FRACTIONS.length;
        if(incidentK==null||incidentK.length!=3)throw new IllegalArgumentException("Actual incident wavevector required");
        double actualKLength=Math.sqrt(sq(incidentK[0],incidentK[1])+incidentK[2]*incidentK[2]);
        if(!Double.isFinite(actualKLength)||Math.abs(actualKLength-k0)>.1||incidentK[2]>=0)throw new IllegalArgumentException("Actual downward vacuum wavevector magnitude required");
        // Surface quadrature can multiply all three constant k components by
        // the same tiny area error. Normalize that measured direction to the
        // exact specified vacuum k0; no wave number is inferred from field data.
        double kx=incidentK[0]*k0/actualKLength,ky=incidentK[1]*k0/actualKLength,kpar=Math.hypot(kx,ky),ct=-incidentK[2]/actualKLength,kz=k0*ct;
        if(!Double.isFinite(kz)||ct<=0||ct>1)throw new IllegalArgumentException("Actual downward nongrazing air wavevector required");
        double px=kpar>.1?kx/kpar:1,py=kpar>.1?ky/kpar:0;
        String dephase="exp(i*(("+kx+"[1/m])*x+("+ky+"[1/m])*y))";
        double[] z=new double[n];double[][]ex=new double[n][2],ey=new double[n][2];boolean areasOK=true;List<Map<String,Object>>planes=new ArrayList<>();
        for(int j=0;j<n;j++){
            z[j]=zland+POWER_FIT_CAP_FRACTIONS[j]*cap;String data="powerFitPlane"+j;
            evaluationCheckpoint("EVALUATING_INTERNAL_E_FIT_PLANE",obj("plane_index",j,"cap_fraction",POWER_FIT_CAP_FRACTIONS[j],"z_m",z[j]));
            m.result().dataset().create(data,"CutPlane");
            try{
                m.result().dataset(data).set("data","dset1");m.result().dataset(data).set("planetype","quick");m.result().dataset(data).set("quickplane","xy");m.result().dataset(data).set("quickz",z[j]+"[m]");
                double a=complexCutIntegral(m,data,"1")[0];areasOK&=close(a,area,GEOMETRY_MEASURE_REL);
                double[]x=complexCutIntegral(m,data,"emw.Ex*("+dephase+")"),y=complexCutIntegral(m,data,"emw.Ey*("+dephase+")");
                ex[j]=new double[]{(px*x[0]+py*y[0])/a,(px*x[1]+py*y[1])/a};ey[j]=new double[]{(-py*x[0]+px*y[0])/a,(-py*x[1]+px*y[1])/a};
                Map<String,Object>r=obj("plane_index",j,"cap_fraction",POWER_FIT_CAP_FRACTIONS[j],"z_m",z[j],"area_m2",a,"relative_area_difference",Math.abs(a-area)/area,"avg_Bloch_dephased_Ep_t_re_im",ex[j],"avg_Bloch_dephased_Es_re_im",ey[j],"avg_Ex_re_im",new double[]{x[0]/a,x[1]/a},"avg_Ey_re_im",new double[]{y[0]/a,y[1]/a},"dephased",true,"phase_convention","exp(+i*omega*t); multiply exp(+i*k_parallel dot r_parallel)","Bloch_dephase_expression",dephase);planes.add(r);evaluationCheckpoint("INTERNAL_E_FIT_PLANE_EVALUATED",r);
            }finally{m.result().dataset().remove(data);}
        }
        double[]fx=fitComplexPlaneWaves(z,ex,kz),fy=fitComplexPlaneWaves(z,ey,kz);
        double pin=tangentialPlanePower(area,z0,ct,new double[]{fx[0],fx[1]},new double[]{fy[0],fy[1]}),pref=tangentialPlanePower(area,z0,ct,new double[]{fx[2],fx[3]},new double[]{fy[2],fy[3]});
        double residual=Math.sqrt((fx[4]+fy[4])/Math.max(fx[5]+fy[5],1e-300)),condition=Math.max(fx[6],fy[6]),maxCondition=condition,maxPinError=0;List<Map<String,Object>>loo=new ArrayList<>();
        for(int held=0;held<n;held++){
            double[]zs=new double[n-1];double[][]xs=new double[n-1][2],ys=new double[n-1][2];int p=0;
            for(int j=0;j<n;j++)if(j!=held){zs[p]=z[j];xs[p]=ex[j];ys[p]=ey[j];p++;}
            double[]xx=fitComplexPlaneWaves(zs,xs,kz),yy=fitComplexPlaneWaves(zs,ys,kz),predP=fittedValue(xx,kz,z[held]),predS=fittedValue(yy,kz,z[held]);
            double lp=tangentialPlanePower(area,z0,ct,new double[]{xx[0],xx[1]},new double[]{yy[0],yy[1]}),err=Math.abs(lp-1),cond=Math.max(xx[6],yy[6]);maxPinError=Math.max(maxPinError,err);maxCondition=Math.max(maxCondition,cond);
            double heldError2=sq(predP[0]-ex[held][0],predP[1]-ex[held][1])+sq(predS[0]-ey[held][0],predS[1]-ey[held][1]);
            loo.add(obj("held_out_plane_index",held,"Pin_W",lp,"Pin_relative_error_vs_specified_1W",err,"design_condition_2",cond,
                "held_out_relative_E_residual",Math.sqrt(heldError2/Math.max(sq(ex[held][0],ex[held][1])+sq(ey[held][0],ey[held][1]),1e-300)),
                "held_out_E_error_over_full_fit_incident_amplitude",Math.sqrt(heldError2/Math.max(sq(fx[0],fx[1])+sq(fy[0],fy[1]),1e-300))));
        }
        Map<String,Object>r=obj("proof_schema",PROOF_SCHEMA,"method",INCIDENT_POWER_METHOD,"scope","Uniform air cap; all five predetermined planes, actual Bloch transverse dephasing, tangential p/s electric components, equal weights, fixed vacuum k0 and kz",
            "actual_k_in_rad_per_m",incidentK,"direction_normalized_to_exact_vacuum_k0",true,"projection_k_in_rad_per_m",new double[]{kx,ky,-kz},"actual_k_magnitude_rad_per_m",actualKLength,"k0_fixed_rad_per_m",k0,"kz_fixed_rad_per_m",kz,"cos_theta",ct,"power_admittance_p_over_vacuum",1/ct,"power_admittance_s_over_vacuum",ct,"tangential_p_axis_xy",new double[]{px,py},"Bloch_dephase_expression",dephase,"cap_fractions",POWER_FIT_CAP_FRACTIONS,"plane_data",planes,"a_p_t_re_im",new double[]{fx[0],fx[1]},"a_s_re_im",new double[]{fy[0],fy[1]},"b_p_t_re_im",new double[]{fx[2],fx[3]},"b_s_re_im",new double[]{fy[2],fy[3]},
            "Pin_W",pin,"Pref_W",pref,"design_condition_2",condition,"maximum_full_or_LOO_design_condition_2",maxCondition,"relative_L2_E_fit_residual",residual,"LOO_max_Pin_relative_error_vs_specified_1W",maxPinError,"leave_one_out",loo,
            "all_cut_areas_match_geometry",areasOK,"area_relative_tolerance",GEOMETRY_MEASURE_REL,"condition_limit",POWER_FIT_CONDITION_MAX,"relative_L2_residual_limit",POWER_FIT_RELATIVE_L2_MAX,"LOO_Pin_error_limit",POWER_FIT_LOO_PIN_ERROR_MAX,
            "condition_passed",maxCondition<=POWER_FIT_CONDITION_MAX,"residual_passed",residual<=POWER_FIT_RELATIVE_L2_MAX,"LOO_Pin_passed",maxPinError<=POWER_FIT_LOO_PIN_ERROR_MAX,"held_out_residual_is_diagnostic_only",true,"R_A_flux_normalization","specified_1W_not_fitted_Pin");
        evaluationCheckpoint("INTERNAL_E_POWER_FIT_EVALUATED",r);return r;
    }
    private static String scatteringVariable(Model m,String pn) {
        String compact="emw.S"+pn+"1",separated="emw.S"+pn+"_1";
        // Multi-digit names differ among interfaces/releases. Ask the actual
        // compiled RF equations; never treat a missing channel as zero power.
        boolean a=global(m,"isdefined("+compact+")")>.5,b=global(m,"isdefined("+separated+")")>.5;
        if(!a&&!b)throw new IllegalStateException("No defined S variable for physical port "+pn);
        if(a&&b&&Math.abs(global(m,"abs("+compact+"-"+separated+")"))>1e-10)throw new IllegalStateException("Ambiguous unequal S aliases for "+pn);
        return b?separated:compact;
    }
    private static Map<String,Map<String,String>> orderDescriptors(PhysicsFeature port) {
        // addDiffractionOrders can delete and recreate its children. Snapshot
        // immutable scalar metadata now; never retain a native feature handle
        // across regeneration, mesh creation, equation compilation or solving.
        Map<String,Map<String,String>> out=new LinkedHashMap<>();
        for(String t:port.feature().tags()){
            PhysicsFeature f=port.feature(t);String type=f.getType();
            if(!type.equals("DiffractionOrder")&&!type.equals("OrthogonalPolarization"))continue;
            boolean orth=type.equals("OrthogonalPolarization");
            String mm=orth?"0":f.getString("mOrder"),nn=orth?"0":f.getString("nOrder"),component=orth?"orthogonal":f.getString("components");
            String key=mm+","+nn+","+component;
            Map<String,String> d=new LinkedHashMap<>();d.put("feature_tag",t);d.put("type",type);
            d.put("port_name",f.getString("PortName"));d.put("m",mm);d.put("n",nn);d.put("components",component);
            if(out.put(key,Collections.unmodifiableMap(d))!=null)throw new IllegalStateException("Duplicate port order "+key);
        }
        return Collections.unmodifiableMap(out);
    }
    public static List<Map<String,Object>> expectedOrders(double lambda,double incidentKx,double incidentKy) {
        List<Map<String,Object>> r=new ArrayList<>();double k=2*Math.PI/(lambda*UM),bx=2*Math.PI/PITCH,by=bx/Math.sqrt(3);
        int bound=(int)Math.ceil(2*(k+Math.hypot(incidentKx,incidentKy))/Math.min(bx,by))+2;
        if(bound>100)throw new IllegalArgumentException("Outside authorized finite reciprocal-lattice enumeration");
        for(int a=-bound;a<=bound;a++)for(int b=-bound;b<=bound;b++){double kx=incidentKx+a*bx,ky=incidentKy+(-a+2*b)*by,q=(kx*kx+ky*ky)/(k*k);if(q<1-1e-12)r.add(obj("m",a,"n",b,"kx_rad_per_m",kx,"ky_rad_per_m",ky,"kz_over_k0",Math.sqrt(1-q)));}
        return r;
    }
    private static boolean validPhysicalFeaturePairs(Map<String,Map<String,String>> fs,int spatialOrders) {
        Map<String,Set<String>> pairs=new LinkedHashMap<>();
        for(Map<String,String> d:fs.values())if(d.get("type").equals("DiffractionOrder")){
            String q=d.get("m")+","+d.get("n");
            pairs.computeIfAbsent(q,k->new HashSet<>()).add(d.get("components"));
        }
        if(pairs.size()!=spatialOrders-1||pairs.containsKey("0,0"))return false;
        for(Set<String>s:pairs.values())if(s.size()!=2||!s.contains("inplane")||!s.contains("outofplane"))return false;
        return fs.containsKey("0,0,orthogonal");
    }
    private static int expectedIndex(List<Map<String,Object>> expected,double kx,double ky) {
        for(int j=0;j<expected.size();j++){Map<String,Object> e=expected.get(j);if(Math.hypot(kx-((Number)e.get("kx_rad_per_m")).doubleValue(),ky-((Number)e.get("ky_rad_per_m")).doubleValue())<.1)return j;}
        return -1;
    }

    /* CODEX_FINITE_COMPLETE_PORT_BASIS_REPAIR_BEGIN */
    private static String canonicalOrderKey(int m,int n){return m+","+n;}
    private static Set<String> independentOrderKeys(List<Map<String,Object>> orders){
        Set<String>s=new LinkedHashSet<>();for(Map<String,Object>r:orders)if(!s.add(canonicalOrderKey(((Number)r.get("m")).intValue(),((Number)r.get("n")).intValue())))throw new IllegalStateException("Duplicate independent reciprocal order");return s;
    }
    private static String nativeDescriptorCanonicalKey(Map<String,String>d){
        return canonicalOrderKey(-Integer.parseInt(d.get("m")),-Integer.parseInt(d.get("n")));
    }
    private static boolean expectedPairsExactly(Map<String,Map<String,String>> fs,Set<String> expected){
        if(!expected.contains("0,0")||fs.size()+1!=2*expected.size()||!fs.containsKey("0,0,orthogonal"))return false;
        Map<String,Set<String>>pairs=new LinkedHashMap<>();Set<String>names=new HashSet<>();names.add("1");
        for(Map<String,String>d:fs.values()){
            if(!names.add(d.get("port_name")))return false;
            String key=nativeDescriptorCanonicalKey(d);if(!expected.contains(key))return false;
            if("OrthogonalPolarization".equals(d.get("type"))){if(!"0,0".equals(key)||!"orthogonal".equals(d.get("components")))return false;continue;}
            if(!"DiffractionOrder".equals(d.get("type"))||"0,0".equals(key))return false;
            if(!pairs.computeIfAbsent(key,k->new HashSet<>()).add(d.get("components")))return false;
        }
        if(pairs.size()!=expected.size()-1)return false;
        for(String key:expected)if(!"0,0".equals(key)){Set<String>s=pairs.get(key);if(s==null||s.size()!=2||!s.contains("inplane")||!s.contains("outofplane"))return false;}
        return true;
    }
    private static Map<String,Object> presolvePortWavevectorProof(Model m,int[] top,double topArea,double lambda,double theta,double phi,double[] requestedK,double[] actualIncidentK,
            List<Map<String,Object>>expected,List<Map<String,Object>>padding,Map<String,Map<String,String>>channels,Set<String>physicalKeys,List<Map<String,Object>>removed){
        progress("PRESOLVE_VERIFYING_COMPLETE_PORT_BASIS");
        double k0=2*Math.PI/(lambda*UM);boolean sourceOK=true;for(int j=0;j<3;j++)sourceOK&=Double.isFinite(actualIncidentK[j])&&Math.abs(actualIncidentK[j]-requestedK[j])<.1;
        sourceOK&=Math.abs(Math.sqrt(actualIncidentK[0]*actualIncidentK[0]+actualIncidentK[1]*actualIncidentK[1]+actualIncidentK[2]*actualIncidentK[2])-k0)<.1&&actualIncidentK[2]<0;
        Set<String>physicalSet=independentOrderKeys(expected),paddingSet=independentOrderKeys(padding);boolean physicalExact=true,paddingExact=expectedPairsExactly(channels,paddingSet),allOK=sourceOK&&paddingExact,namesUnique=true;
        Map<String,Integer>physicalCount=new LinkedHashMap<>(),paddingCount=new LinkedHashMap<>();Set<String>names=new HashSet<>();List<Map<String,Object>>rows=new ArrayList<>();
        List<Map<String,String>>descriptors=new ArrayList<>();Map<String,String>zero=new LinkedHashMap<>();zero.put("port_name","1");zero.put("m","0");zero.put("n","0");zero.put("components","co");zero.put("key","0,0,co");descriptors.add(zero);
        for(Map.Entry<String,Map<String,String>>e:channels.entrySet()){Map<String,String>d=new LinkedHashMap<>(e.getValue());d.put("key",e.getKey());descriptors.add(d);}
        int[]spatialCounts=new int[expected.size()];
        for(Map<String,String>d:descriptors){
            String pn=d.get("port_name"),key=d.get("key");int nativeM=Integer.parseInt(d.get("m")),nativeN=Integer.parseInt(d.get("n")),cm=-nativeM,cn=-nativeN;String ck=canonicalOrderKey(cm,cn);boolean prop=physicalSet.contains(ck);
            if(!names.add(pn))namesUnique=false;
            double kx=avg(m,top,topArea,"real(emw.kModex_"+pn+")"),ky=avg(m,top,topArea,"real(emw.kModey_"+pn+")"),kz=avg(m,top,topArea,"real(emw.kModez_"+pn+")"),ikz=avg(m,top,topArea,"imag(emw.kModez_"+pn+")");
            double ex=requestedK[0]+cm*2*Math.PI/PITCH,ey=requestedK[1]+(-cm+2*cn)*2*Math.PI/PITCH/Math.sqrt(3),q=(kx*kx+ky*ky)/(k0*k0);boolean actualProp=q<1-1e-10;
            int expectedId=expectedIndex(expected,kx,ky);boolean finite=Double.isFinite(kx)&&Double.isFinite(ky)&&Double.isFinite(kz)&&Double.isFinite(ikz);
            boolean waveOK=finite&&Math.abs(kx-ex)<.1&&Math.abs(ky-ey)<.1&&actualProp==prop&&paddingSet.contains(ck)&&Math.abs((kx*kx+ky*ky+kz*kz-ikz*ikz)/(k0*k0)-1)<1e-7;
            if(prop){waveOK&=expectedId>=0&&Math.abs(ikz)<.1;if(expectedId>=0)spatialCounts[expectedId]++;physicalCount.put(ck,physicalCount.getOrDefault(ck,0)+1);}
            // Compare axial decay with the exact independent reference, not a
            // near-cutoff square root of area-quadrature-noisy measured kx/ky.
            // Raw measured components and dispersion retain their .1/1e-7 gates.
            else waveOK&=expectedId<0&&q>1+1e-10&&Math.abs(kz)<.1&&Math.abs(Math.abs(ikz)-k0*Math.sqrt((ex*ex+ey*ey)/(k0*k0)-1))<.1;
            paddingCount.put(ck,paddingCount.getOrDefault(ck,0)+1);allOK&=waveOK;
            rows.add(obj("port_name",pn,"key",key,"native_m",nativeM,"native_n",nativeN,"canonical_m",cm,"canonical_n",cn,"components",d.get("components"),"propagating",prop,"actual_wavevector_propagating",actualProp,"independent_expected_spatial_index",expectedId,"expected_kx_rad_per_m",ex,"expected_ky_rad_per_m",ey,"kx_rad_per_m",kx,"ky_rad_per_m",ky,"real_kz_rad_per_m",kz,"imag_kz_rad_per_m",ikz,"wavevector_verified",waveOK));
        }
        physicalExact=physicalCount.keySet().equals(physicalSet);for(int n:physicalCount.values())physicalExact&=n==2;for(int n:spatialCounts)physicalExact&=n==2;
        paddingExact&=paddingCount.keySet().equals(paddingSet);for(int n:paddingCount.values())paddingExact&=n==2;
        int physicalChannels=physicalKeys.size()+1,total=channels.size()+1;allOK&=physicalExact&&paddingExact&&namesUnique&&physicalChannels==2*expected.size()&&total==2*padding.size();
        Map<String,Object>proof=obj("case_id",caseId,"proof_schema",PROOF_SCHEMA,"state",allOK?"PASS_PRESOLVE_COMPLETE_PORT_BASIS":"FAIL_PRESOLVE_COMPLETE_PORT_BASIS","solver_executed",false,"physical_frequency_restored",true,"source_frequency_expression",m.study("std1").feature("freq").getString("plist"),"generator_superset_frequency_ratio",1.20,"retained_padding_frequency_ratio",1.10,"wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"requested_k_in_rad_per_m",requestedK,"actual_incident_k_rad_per_m",actualIncidentK,"incident_wavevector_matches_requested",sourceOK,"actual_wavevectors_match_independent_lattice",allOK,"physical_spatial_orders",expected.size(),"physical_polarization_channels",physicalChannels,"total_port_channels",total,"evanescent_channels",total-physicalChannels,"actual_independent_spatial_channel_counts",spatialCounts,"physical_index_set_exact",physicalExact,"padding_index_set_exact",paddingExact,"port_names_unique",namesUnique,"vector_absolute_tolerance_rad_per_m",.1,"dispersion_relative_tolerance",1e-7,"native_to_canonical_index_mapping","canonical(m,n)=native(-m,-n); verified independently by every actual global wavevector","independent_expected_lattice_orders",expected,"independent_retained_padding_orders",padding,"channels",rows,"removed_superset_channels",removed,"native_field_solver_run_before_proof",false,"updated_utc",Instant.now().toString());
        emit("CODEX_PRESOLVE_PORT_BASIS_JSON",proof);if(!allOK||!"f0".equals(proof.get("source_frequency_expression")))throw new IllegalStateException("Pre-solve complete physical/padding wavevector admission failed");return proof;
    }
    /* CODEX_FINITE_COMPLETE_PORT_BASIS_REPAIR_END */

    private static String triKey(double[][] xyz,int[][] tri,int idx,double dx,double dy) {
        String[] p=new String[3];for(int j=0;j<3;j++){int n=tri[j][idx];p[j]=Math.round((xyz[0][n]+dx)/1e-12)+":"+Math.round((xyz[1][n]+dy)/1e-12)+":"+Math.round(xyz[2][n]/1e-12);}
        Arrays.sort(p);return p[0]+"/"+p[1]+"/"+p[2];
    }
    private static SortedMap<String,Integer> faceTriangles(double[][] xyz,int[][] tri,int[] entities,int face,double dx,double dy) {SortedMap<String,Integer>s=new TreeMap<>();for(int i=0;i<entities.length;i++)if(entities[i]==face){String k=triKey(xyz,tri,i,dx,dy);s.put(k,s.getOrDefault(k,0)+1);}return s;}
    private static String segmentMask(double[] a,double[] b) {
        String[] r=new String[3];String[] axes={"x","y","z"};double[] d=new double[3];double d2=0;
        StringJoiner dot=new StringJoiner("+");
        for(int j=0;j<3;j++){r[j]="("+axes[j]+"-("+a[j]+"[m]))";d[j]=b[j]-a[j];d2+=d[j]*d[j];dot.add(r[j]+"*("+d[j]+"[m])");}
        String u="min(1,max(0,("+dot+")/("+d2+"[m^2])))";
        StringJoiner dist2=new StringJoiner("+");for(int j=0;j<3;j++)dist2.add("("+r[j]+"-("+u+")*("+d[j]+"[m]))^2");
        return "("+dist2+"<edgeStrip^2)";
    }
    private static String edgeMask(boolean flat,double HEIGHT) {
        if(flat)return "0";
        // True Euclidean distance to the 18 finite sharp edge segments, not
        // infinite side lines (which would overcount land near hex corners).
        double[][] b=new double[6][3],t=new double[6][3];
        for(int j=0;j<6;j++){double q=Math.PI/6+j*Math.PI/3;b[j]=new double[]{ABOTTOM/Math.cos(Math.PI/6)*Math.cos(q),ABOTTOM/Math.cos(Math.PI/6)*Math.sin(q),0};t[j]=new double[]{AOPEN/Math.cos(Math.PI/6)*Math.cos(q),AOPEN/Math.cos(Math.PI/6)*Math.sin(q),HEIGHT};}
        StringJoiner union=new StringJoiner("||","(",")");
        for(int j=0;j<6;j++){union.add(segmentMask(b[j],b[(j+1)%6]));union.add(segmentMask(t[j],t[(j+1)%6]));union.add(segmentMask(b[j],t[j]));}
        return "if("+union+",1,0)";
    }

    private static void digestInt(MessageDigest d,int n){for(int j=3;j>=0;j--)d.update((byte)(n>>>(j*8)));}
    private static void digestLong(MessageDigest d,long n){for(int j=7;j>=0;j--)d.update((byte)(n>>>(j*8)));}
    private static void digestInts(MessageDigest d,int[]a){digestInt(d,a.length);for(int n:a)digestInt(d,n);}
    private static void digestMatrix(MessageDigest d,int[][]a){digestInt(d,a.length);for(int[]row:a)digestInts(d,row);}
    /** Pure in-memory exact mesh identity: includes ordering, not merely counts.
     * SHA256 over version tag, big-endian dimensions, IEEE754 coordinate bits,
     * tetrahedral and triangular connectivity, domain and boundary entity IDs.
     */
    private static String meshFingerprint(MeshSequence ms)throws Exception{
        MessageDigest d=MessageDigest.getInstance("SHA-256");d.update(new byte[]{67,80,77,69,83,72,1});
        double[][]v=ms.getVertex();digestInt(d,v.length);for(double[]row:v){digestInt(d,row.length);for(double x:row){if(!Double.isFinite(x))throw new IllegalStateException("Nonfinite mesh coordinate");digestLong(d,Double.doubleToLongBits(x));}}
        digestMatrix(d,ms.getElem("tet"));digestInts(d,ms.getElemEntity("tet"));digestMatrix(d,ms.getElem("tri"));digestInts(d,ms.getElemEntity("tri"));
        StringBuilder s=new StringBuilder();for(byte b:d.digest()){s.append("0123456789abcdef".charAt((b>>>4)&15));s.append("0123456789abcdef".charAt(b&15));}return s.toString();
    }

    private static double volumeIntegral(Model m,int[] domains,String expr) {
        String t="ev"+(++evalId);lastNumericalTag=t;lastNumericalKind="IntVolume";lastExpression=expr;m.result().numerical().create(t,"IntVolume");
        try{NumericalFeature q=m.result().numerical(t);q.selection().set(domains);q.set("data","dset1");q.set("expr",expr);q.set("intorderactive","on");q.set("intorder",8);double r=q.getReal()[0][0];if(!Double.isFinite(r))throw new IllegalStateException("Nonfinite volume loss");return r;}
        finally{m.result().numerical().remove(t);}
    }

    /* CODEX_EXACT_OOC_QUEUE_CONFIGURATION_BEGIN */
    private static int authorizedIndex(String id) {
        String[] ids={"ML6_A3H100_t2_l10p5_th0_phi0_s","ML7_A3H100_t2_l10p5_th0_phi0_p","ML8_A3H100_t2_l10p5_th30_phi30_s","ML9_A3H100_t2_l10p5_th30_phi30_p"};for(int i=0;i<ids.length;i++)if(ids[i].equals(id))return i;
        throw new IllegalArgumentException("Only four closed ML remaining H100 memory-recovery IDs exist");
    }
    private static String authorizedPolarization(String id){String[]v={"s","p","s","p"};return v[authorizedIndex(id)];}
    private static int authorizedFilmUm(String id){int[]v={2,2,2,2};return v[authorizedIndex(id)];}
    private static double authorizedLambda(String id){double[]v={10.5,10.5,10.5,10.5};return v[authorizedIndex(id)];}
    private static double authorizedTheta(String id){double[]v={0,0,30,30};return v[authorizedIndex(id)];}
    private static double authorizedPhi(String id){double[]v={0,0,30,30};return v[authorizedIndex(id)];}
    private static double authorizedHeightUm(String id){double[]v={100,100,100,100};return v[authorizedIndex(id)];}
    private static String authorizedMesh(String id){String[]v={"m4_air28","m4_air28","m4_air28","m4_air28"};return v[authorizedIndex(id)];}
    private static int authorizedPaddingChannels(String id){int[]v={398,398,380,380};return v[authorizedIndex(id)];}
    private static double[] authorizedOptical(String id){double[][]v={{2.2447212969283283,0.06982232081911283,4.308219763636363,20.63042690909091},{2.2447212969283283,0.06982232081911283,4.308219763636363,20.63042690909091},{2.2447212969283283,0.06982232081911283,4.308219763636363,20.63042690909091},{2.2447212969283283,0.06982232081911283,4.308219763636363,20.63042690909091}};return v[authorizedIndex(id)];}
    /** Incoming and outgoing true global p; both are transverse to their k. */
    public static double[] globalP(double theta,double phi,boolean outgoing){double t=Math.toRadians(theta),p=Math.toRadians(phi);return new double[]{Math.cos(t)*Math.cos(p),Math.cos(t)*Math.sin(p),(outgoing?-1:1)*Math.sin(t)};}
    /** True global s is perpendicular to the global azimuth plane; never merely relabel p. */
    public static double[] globalE(String pol,double theta,double phi,boolean outgoing){if("p".equals(pol))return globalP(theta,phi,outgoing);if("s".equals(pol)){double p=Math.toRadians(phi);return new double[]{-Math.sin(p),Math.cos(p),0};}throw new IllegalArgumentException("Only explicit true s/p polarization");}
    /** Exact direct factorization with disk-backed LU blocks, not BLR or a
     * process-memory guarantee. Never changes the active fc1 -> d1 selection.
     * apply=true is called before mesh_only.mph is saved; apply=false performs
     * a second strong readback immediately before the physical solver call.
     */
    private static Map<String,Object> exactOutOfCore(Model m,String solver,boolean apply) {
        authorizedPolarization(caseId);
        SolverFeature stationary=m.sol(solver).feature("s1");
        if(!"Stationary".equals(stationary.getType()))throw new IllegalStateException("Expected unchanged Stationary s1");
        int stationaryCount=0;for(String t:m.sol(solver).feature().tags())if("Stationary".equals(m.sol(solver).feature(t).getType()))stationaryCount++;
        if(stationaryCount!=1)throw new IllegalStateException("Unique Stationary solver required");
        Map<String,Object> activeBefore=new LinkedHashMap<>();List<String> directTags=new ArrayList<>();
        for(String tag:stationary.feature().tags()){
            SolverFeature f=stationary.feature(tag);
            if("FullyCoupled".equals(f.getType()))activeBefore.put(tag,f.getString("linsolver"));
            if("Direct".equals(f.getType()))directTags.add(tag);
        }
        if(activeBefore.size()!=1||!"d1".equals(activeBefore.get("fc1"))||directTags.size()!=2||!directTags.contains("dDef")||!directTags.contains("d1"))throw new IllegalStateException("Expected original fc1 -> d1 and exactly dDef/d1 Direct features");
        List<Map<String,Object>> changes=new ArrayList<>();
        for(String tag:directTags){
            SolverFeature d=stationary.feature(tag);String expectedLinear="dDef".equals(tag)?"mumps":"pardiso",linear=d.getString("linsolver");
            if(!expectedLinear.equals(linear))throw new IllegalStateException("Original Direct algorithm changed: "+tag);
            Map<String,Object> before=obj("linsolver",linear,"ooc",d.getString("ooc"),"incore",d.getString("incore"),"oocmemory_MB",d.getDouble("oocmemory"),"mumpsblr",d.getString("mumpsblr"));
            if(apply){d.set("ooc","on");d.set("incore","manual");d.set("oocmemory",24576.0);d.set("mumpsblr","off");}
            Map<String,Object> after=obj("linsolver",d.getString("linsolver"),"ooc",d.getString("ooc"),"incore",d.getString("incore"),"oocmemory_MB",d.getDouble("oocmemory"),"mumpsblr",d.getString("mumpsblr"));
            if(!linear.equals(d.getString("linsolver"))||!"on".equals(d.getString("ooc"))||!"manual".equals(d.getString("incore"))||Double.compare(d.getDouble("oocmemory"),24576.0)!=0||!"off".equals(d.getString("mumpsblr")))throw new IllegalStateException("Exact OOC24GiB readback failed: "+tag);
            changes.add(obj("tag",tag,"before",before,"after",after,"readback_verified",true));
        }
        Map<String,Object> activeAfter=new LinkedHashMap<>();for(String tag:stationary.feature().tags())if("FullyCoupled".equals(stationary.feature(tag).getType()))activeAfter.put(tag,stationary.feature(tag).getString("linsolver"));
        if(!activeBefore.equals(activeAfter))throw new IllegalStateException("Active linear solver selection changed");
        return obj("case_id",caseId,"state","VERIFIED_OOC_CONFIGURATION_BEFORE_PHYSICS","method","EXACT_DIRECT_OUT_OF_CORE_NO_BLR","solver_tag",solver,"stationary_tag","s1","active_selection_before",activeBefore,"active_selection_after",activeAfter,"active_selection_preserved",true,"direct_feature_count",changes.size(),"direct_changes",changes,"readback_verified",true,"oocmemory_is_not_process_RSS_cap",true,"physical_thresholds_changed",false);
    }
    /* CODEX_EXACT_OOC_QUEUE_CONFIGURATION_END */

    public static void run(String absoluteCaseDir,String id,double nfilm,double kfilm,String mesh,boolean meshOnly) throws Exception {
        root=absoluteCaseDir.replace('/','\\').replace("\\\\","\\");caseId=id;started=Instant.now().toString();evalId=0;
        lastPhase="INITIALIZING";lastExpression=null;lastNumericalTag=null;lastNumericalKind=null;solverCompleted=false;solvedArtifactSaved=false;
        if(!root.startsWith("D:\\meep_sim_evidence\\")||root.contains("..")||!id.matches("[A-Za-z0-9_]+"))throw new IllegalArgumentException("Owned explicit directory/id required");
        if(!Double.isFinite(nfilm)||!Double.isFinite(kfilm)||nfilm<=0||kfilm<=0||!authorizedMesh(id).equals(mesh))throw new IllegalArgumentException("Passive film and exact id-bound mesh required");
        final int filmUm=authorizedFilmUm(id);final double FILM_THICKNESS=filmUm*UM;
        final double HEIGHT=authorizedHeightUm(id)*UM;
        // Closed five-case release: no8um or arbitrary source/mesh fallback.
        boolean flat=false,pec=false;String geometry="A3_FLOOR_FILL_"+filmUm+"UM_ONLY",surface="TI_IBC_WITH_SIO2_FLOOR_FILL",pol=authorizedPolarization(id);
        double lambda=authorizedLambda(id),theta=authorizedTheta(id),phi=authorizedPhi(id),cap=10.5*UM,zland=HEIGHT,ztop=HEIGHT+cap;
        double[] nk=authorizedOptical(id);double nti=nk[2],kti=nk[3];
        double expectedNFilm=nk[0],expectedKFilm=nk[1];
        if(Math.abs(nfilm-expectedNFilm)>1e-14||Math.abs(kfilm-expectedKFilm)>1e-14)throw new IllegalArgumentException("Film source n/k does not match the authorized wavelength interpolation");
        double[] requestedK=requestedIncidentK(lambda,theta,phi);double tr=Math.toRadians(theta);
        double[] sourcePlaneZUm=new double[POWER_FIT_CAP_FRACTIONS.length];for(int j=0;j<sourcePlaneZUm.length;j++)sourcePlaneZUm[j]=(HEIGHT+POWER_FIT_CAP_FRACTIONS[j]*cap)/UM;
        double[] incomingE=globalE(pol,theta,phi,false),outgoingE=globalE(pol,theta,phi,true);
        Model m=null;
        try {
            progress("BUILDING_GEOMETRY");m=ModelUtil.create("CavityFloorFilm");m.label("A3 explicit"+filmUm+"um SiO2 floor-fill controlled pilot; NOT all-wall coating");
            m.param().set("lambda0",lambda+"[um]");m.param().set("f0","c_const/lambda0");m.param().set("Pin0","1[W]");
            m.param().set("nti",Double.toString(nti));m.param().set("kti",Double.toString(kti));m.param().set("Zs","Z0_const/(nti-i*kti)");
            m.param().set("Hc",authorizedHeightUm(id)+"[um]");m.param().set("ab","7.5[um]");m.param().set("ao","34.77279213156495[um]");m.param().set("edgeStrip","0.5[um]");m.param().set("tol","1e-5[um]");
            if(!(authorizedHeightUm(id)+"[um]").equals(m.param().get("Hc")))throw new IllegalStateException("Case-bound Hc parameter readback changed");
            m.param().set("R_ibc",pol.equals("s")?"abs((Zs*"+Math.cos(tr)+"-Z0_const)/(Zs*"+Math.cos(tr)+"+Z0_const))^2":"abs((Zs-Z0_const*"+Math.cos(tr)+")/(Zs+Z0_const*"+Math.cos(tr)+"))^2");
            m.param().set("nfilm",Double.toString(nfilm));m.param().set("kfilm",Double.toString(kfilm));m.param().set("eps_f","(nfilm-i*kfilm)^2");
            String nfReadback=m.param().get("nfilm"),kfReadback=m.param().get("kfilm"),epsReadback=m.param().get("eps_f");
            String ntiReadback=m.param().get("nti"),ktiReadback=m.param().get("kti"),zsReadback=m.param().get("Zs");
            boolean titaniumParameterBinding=Double.toString(nti).equals(ntiReadback)&&Double.toString(kti).equals(ktiReadback)&&"Z0_const/(nti-i*kti)".equals(zsReadback);
            boolean filmParameterBinding=Double.toString(nfilm).equals(nfReadback)&&Double.toString(kfilm).equals(kfReadback)&&"(nfilm-i*kfilm)^2".equals(epsReadback)&&titaniumParameterBinding;
            emit("CODEX_MATERIAL_CONFIGURATION_JSON",obj("proof_schema",PROOF_SCHEMA,"wavelength_um",lambda,"nfilm",nfilm,"kfilm",kfilm,"nfilm_readback",nfReadback,"kfilm_readback",kfReadback,"eps_f_readback",epsReadback,"eps_f_real",nfilm*nfilm-kfilm*kfilm,"eps_f_imag",-2*nfilm*kfilm,"nti",nti,"kti",kti,"nti_readback",ntiReadback,"kti_readback",ktiReadback,"Zs_readback",zsReadback,"titanium_parameter_binding_verified",titaniumParameterBinding,"film_parameter_binding_verified",filmParameterBinding,"film_model","explicit_lossy_volume_not_layered_IBC"));
            if(!filmParameterBinding)throw new IllegalStateException("Film n/k/epsilon parameter binding mismatch");
            m.component().create("comp1",true);m.component("comp1").geom().create("geom1",3);GeomSequence g=m.component("comp1").geom("geom1");g.lengthUnit("m");
            double fillApothem=ABOTTOM+(AOPEN-ABOTTOM)*FILM_THICKNESS/HEIGHT;
            hex(g,"wpFilm","film",ABOTTOM,0,FILM_THICKNESS,fillApothem/ABOTTOM);
            hex(g,"wpCavity","cavity",fillApothem,FILM_THICKNESS,HEIGHT-FILM_THICKNESS,AOPEN/fillApothem);
            hex(g,"wpCap","cap",APOTHEM,zland,cap,1);
            g.create("airUnion","Union");g.feature("airUnion").selection("input").set(new String[]{"cavity","cap"});g.feature("airUnion").set("intbnd",false);
            g.create("twoMaterials","Union");g.feature("twoMaterials").selection("input").set(new String[]{"film","airUnion"});g.feature("twoMaterials").set("intbnd",true);
            g.run();GeomInfo gi=g;gi.check();double[][] vertices=gi.getVertexCoord();
            int filmDomain=0,airDomain=0;
            for(int d=1;d<=gi.getNDomains();d++){double lo=Double.POSITIVE_INFINITY,hi=Double.NEGATIVE_INFINITY;for(int v:gi.getAdj(3,0,d)){lo=Math.min(lo,vertices[2][v-1]);hi=Math.max(hi,vertices[2][v-1]);}
                if(Math.abs(lo)<TOL&&Math.abs(hi-FILM_THICKNESS)<TOL){if(filmDomain!=0)throw new IllegalStateException("Duplicate film domain");filmDomain=d;}
                else if(Math.abs(lo-FILM_THICKNESS)<TOL&&Math.abs(hi-ztop)<TOL){if(airDomain!=0)throw new IllegalStateException("Duplicate air domain");airDomain=d;}
                else throw new IllegalStateException("Unexpected material-domain extent");}
            if(gi.getNDomains()!=2||filmDomain==0||airDomain==0)throw new IllegalStateException("Exactly lower-frustum film and upper-cavity/cap air required");
            List<Integer> topL=new ArrayList<>(),floorL=new ArrayList<>(),landL=new ArrayList<>(),wallL=new ArrayList<>(),wallFilmL=new ArrayList<>(),wallAirL=new ArrayList<>(),interfaceL=new ArrayList<>();int[] sides=new int[6];
            List<Map<String,Object>> faceAudit=new ArrayList<>();
            for(int f=1;f<=gi.getNFaces();f++){
                int[] vv=gi.getAdj(2,0,f),adj=gi.getAdj(2,3,f);String kind=null;int sideIdx=-1;
                if(allZ(vertices,vv,ztop)){topL.add(f);kind="top";if(adj.length!=1||adj[0]!=airDomain)throw new IllegalStateException("Top is not air-adjacent");}
                else if(allZ(vertices,vv,0)){floorL.add(f);kind="Ti_floor_film_side";if(adj.length!=1||adj[0]!=filmDomain)throw new IllegalStateException("Ti floor not film-adjacent");}
                else if(allZ(vertices,vv,FILM_THICKNESS)){interfaceL.add(f);kind="film_air_internal_continuity";if(adj.length!=2||!((adj[0]==filmDomain&&adj[1]==airDomain)||(adj[1]==filmDomain&&adj[0]==airDomain)))throw new IllegalStateException("Film-air interface adjacency wrong");}
                else if(allZ(vertices,vv,zland)){landL.add(f);kind="Ti_land_bare_air_side";if(adj.length!=1||adj[0]!=airDomain)throw new IllegalStateException("Land not air-adjacent");}
                else {for(int j=0;j<6;j++)if(side(vertices,vv,j*Math.PI/3,zland)){if(sides[j]!=0)throw new IllegalStateException("Duplicate cap side");sides[j]=f;sideIdx=j;kind="periodic_air_side";if(adj.length!=1||adj[0]!=airDomain)throw new IllegalStateException("Periodic side not homogeneous air");break;}
                    if(sideIdx<0){boolean onOriginalWall=false;for(int j=0;j<6;j++){boolean all=vv.length>0;for(int v:vv){double z=vertices[2][v-1],rho=vertices[0][v-1]*Math.cos(j*Math.PI/3)+vertices[1][v-1]*Math.sin(j*Math.PI/3);all&=z>=-TOL&&z<=HEIGHT+TOL&&Math.abs(rho-(ABOTTOM+(AOPEN-ABOTTOM)*z/HEIGHT))<TOL;}onOriginalWall|=all;}
                        if(!onOriginalWall||adj.length!=1)throw new IllegalStateException("Unknown face must never receive Ti IBC");wallL.add(f);
                        if(adj[0]==filmDomain){wallFilmL.add(f);kind="Ti_lower_wall_film_side";}else if(adj[0]==airDomain){wallAirL.add(f);kind="Ti_upper_wall_bare_air_side";}else throw new IllegalStateException("Unknown wall material adjacency");}}
                faceAudit.add(obj("face_id",f,"kind",kind,"side_index",sideIdx,"vertices",vv,"adjacent_domains",adj,"area_m2",area(m,new int[]{f})));
            }
            int[] top=ints(topL),floor=ints(floorL),land=ints(landL),walls=ints(wallL),metal=combine(floor,land,walls),filmTi=combine(floor,ints(wallFilmL)),airTi=combine(land,ints(wallAirL)),filmInterface=ints(interfaceL);
            boolean geometryOK=top.length==1&&floor.length==1&&walls.length==12&&wallFilmL.size()==6&&wallAirL.size()==6&&land.length>0&&filmInterface.length==1;for(int x:sides)geometryOK&=x>0;
            double cellArea=2*Math.sqrt(3)*APOTHEM*APOTHEM,openArea=2*Math.sqrt(3)*AOPEN*AOPEN,bottomArea=2*Math.sqrt(3)*ABOTTOM*ABOTTOM,fillArea=2*Math.sqrt(3)*fillApothem*fillApothem;
            double expectedVolume=cellArea*cap+HEIGHT*(openArea+bottomArea+Math.sqrt(openArea*bottomArea))/3;
            double expectedFilmVolume=FILM_THICKNESS*(bottomArea+fillArea+Math.sqrt(bottomArea*fillArea))/3;
            GeomMeasureFinal gm=g.measureFinal();gm.selection().geom("geom1",3);gm.selection().all();double volume=gm.getVolume();gm.selection().set(new int[]{filmDomain});double filmVolume=gm.getVolume();gm.selection().set(new int[]{airDomain});double airVolume=gm.getVolume();
            double topArea=area(m,top),floorArea=area(m,floor),landArea=area(m,land),expectedWallArea=2*Math.sqrt(3)*(AOPEN+ABOTTOM)*Math.sqrt(HEIGHT*HEIGHT+(AOPEN-ABOTTOM)*(AOPEN-ABOTTOM));
            geometryOK&=close(volume,expectedVolume,GEOMETRY_MEASURE_REL)&&close(filmVolume,expectedFilmVolume,GEOMETRY_MEASURE_REL)&&close(airVolume+filmVolume,volume,GEOMETRY_MEASURE_REL)&&close(topArea,cellArea,GEOMETRY_MEASURE_REL)&&close(floorArea,bottomArea,GEOMETRY_MEASURE_REL)&&close(landArea,cellArea-openArea,GEOMETRY_MEASURE_REL)&&close(area(m,walls),expectedWallArea,GEOMETRY_MEASURE_REL)&&close(area(m,filmInterface),fillArea,GEOMETRY_MEASURE_REL);
            int referencePoint=-1;double[] ref=new double[]{APOTHEM,APOTHEM/Math.sqrt(3),ztop};
            for(int j=0;j<vertices[0].length;j++)if(Math.abs(vertices[0][j]-ref[0])<TOL&&Math.abs(vertices[1][j]-ref[1])<TOL&&Math.abs(vertices[2][j]-ref[2])<TOL)referencePoint=j+1;
            geometryOK&=referencePoint>0;
            emit("CODEX_GEOMETRY_JSON",obj("case_id",caseId,"height_um",HEIGHT/UM,"air_cap_um",cap/UM,"top_plane_z_um",ztop/UM,"fixed_uniform_air_source_plane_z_um",sourcePlaneZUm,"proof_schema",PROOF_SCHEMA,"state","GEOMETRY_CHECKED","geometry_verified",geometryOK,"geometry_measure_relative_tolerance",GEOMETRY_MEASURE_REL,"vertex_coordinate_tolerance_m",TOL,"geometry",geometry,"surface",surface,"units","m","air_volume_m3",airVolume,"film_volume_m3",filmVolume,"expected_film_volume_m3",expectedFilmVolume,"total_computational_volume_m3",volume,"expected_total_volume_m3",expectedVolume,"film_air_interface_area_m2",area(m,filmInterface),"floor_fill_thickness_um",filmUm,"film_top_apothem_um",fillApothem/UM,"film_domain",filmDomain,"air_domain",airDomain,"film_interface_faces",filmInterface,"Ti_film_side_faces",filmTi,"Ti_air_side_faces",airTi,"not_conformal_wall_or_land_coating",true,"top_area_m2",topArea,"floor_area_m2",floorArea,"land_area_m2",landArea,"wall_area_m2",area(m,walls),"faces",faceAudit,"periodic_side_ids",sides,"port_reference_vertex",referencePoint,"port_reference_xyz_m",ref,"global_E_direction",incomingE,"wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"normal_incidence_only",theta==0,"single_point_spectral_angular_exploration",true));
            emit("CODEX_TI_GROUP_GEOMETRY_JSON",obj("Ti_film_side_area_m2",area(m,filmTi),"Ti_air_side_area_m2",area(m,airTi),"Ti_floor_area_m2",floorArea,"Ti_lower_wall_area_m2",area(m,ints(wallFilmL)),"Ti_upper_wall_area_m2",area(m,ints(wallAirL)),"Ti_land_area_m2",landArea,"film_air_interface_area_m2",area(m,filmInterface),"film_air_interface_excluded_from_IBC",true,"film_volume_m3",filmVolume,"air_volume_m3",airVolume));
            if(!geometryOK)throw new IllegalStateException("Actual geometry/areas/port vertex fail approved geometry");
            explicit(m,"top",2,top);explicit(m,"metal",2,metal);explicit(m,"filmDomain",3,new int[]{filmDomain});explicit(m,"airDomain",3,new int[]{airDomain});explicit(m,"filmTi",2,filmTi);explicit(m,"airTi",2,airTi);explicit(m,"filmInterface",2,filmInterface);
            m.component("comp1").material().create("airmat","Common");m.component("comp1").material("airmat").selection().set(new int[]{airDomain});m.component("comp1").material("airmat").propertyGroup("def").set("relpermittivity","1");m.component("comp1").material("airmat").propertyGroup("def").set("relpermeability","1");m.component("comp1").material("airmat").propertyGroup("def").set("electricconductivity","0[S/m]");
            m.component("comp1").material().create("filmmat","Common");m.component("comp1").material("filmmat").selection().set(new int[]{filmDomain});m.component("comp1").material("filmmat").propertyGroup("def").set("relpermittivity","eps_f");m.component("comp1").material("filmmat").propertyGroup("def").set("relpermeability","1");m.component("comp1").material("filmmat").propertyGroup("def").set("electricconductivity","0[S/m]");
            m.component("comp1").physics().create("emw","ElectromagneticWaves","geom1");
            m.component("comp1").physics("emw").prop("ShapeProperty").set("order_electricfield","2");
            String fieldOrder=m.component("comp1").physics("emw").prop("ShapeProperty").getString("order_electricfield");
            if(!fieldOrder.equals("2"))throw new IllegalStateException("Actual electric curl element order is not quadratic");
            m.component("comp1").physics("emw").create("port1","Port",2);PhysicsFeature port=m.component("comp1").physics("emw").feature("port1");
            port.selection().set(top);port.set("PortType","Periodic");port.set("PortName","1");port.set("PortExcitation","on");port.set("Pin","Pin0");
            String[] desiredGlobalE=new String[]{Double.toString(outgoingE[0])+"[V/m]",Double.toString(outgoingE[1])+"[V/m]",Double.toString(outgoingE[2])+"[V/m]"};
            port.set("Polarization","UserDefined");port.set("InputType","E");port.set("E0",desiredGlobalE);port.set("Eampl",desiredGlobalE);port.set("alpha1_inc",theta+"[deg]");port.set("alpha2_inc",(phi+PORT_AZIMUTH_OFFSET_DEG)+"[deg]");
            String[] actualE0=port.getStringArray("E0"),actualEampl=port.getStringArray("Eampl");
            boolean sourceBinding=Arrays.equals(desiredGlobalE,actualE0)&&Arrays.equals(desiredGlobalE,actualEampl);
            Map<String,Object> sourceProof=obj("case_id",caseId,"height_um",HEIGHT/UM,"air_cap_um",cap/UM,"top_plane_z_um",ztop/UM,"fixed_uniform_air_source_plane_z_um",sourcePlaneZUm,"proof_schema",PROOF_SCHEMA,"lambda_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"k_in_requested_rad_per_m",requestedK,"k_out_requested_rad_per_m",new double[]{requestedK[0],requestedK[1],-requestedK[2]},"incoming_E_unit_global",incomingE,"outgoing_E_unit_global",outgoingE,"alpha1_inc_readback",port.getString("alpha1_inc"),"alpha2_inc_readback",port.getString("alpha2_inc"),"port_azimuth_offset_deg",PORT_AZIMUTH_OFFSET_DEG,"global_E_requested",Arrays.asList(desiredGlobalE),"E0_readback",Arrays.asList(actualE0),"Eampl_readback",Arrays.asList(actualEampl),"explicit_global_source_binding_verified",sourceBinding,"mode_direction_leakage_limit",MODE_DIRECTION_LEAKAGE_MAX,"incident_power_method",INCIDENT_POWER_METHOD);
            emit("CODEX_SOURCE_CONFIGURATION_JSON",sourceProof);
            if(!sourceBinding)throw new IllegalStateException("Explicit E0/Eampl global source binding readback mismatch");
            port.create("ref1","PeriodicPortReferencePoint",0);port.feature("ref1").selection().set(new int[]{referencePoint});
            for(int j=0;j<3;j++){String t="pc"+j;m.component("comp1").physics("emw").create(t,"PeriodicCondition",2);PhysicsFeature p=m.component("comp1").physics("emw").feature(t);p.selection().set(new int[]{sides[j],sides[j+3]});p.set("PeriodicType","Floquet");p.set("Floquet_source","FromPeriodicPort");}
            if(pec){m.component("comp1").physics("emw").create("pecTi","PerfectElectricConductor",2);m.component("comp1").physics("emw").feature("pecTi").selection().set(metal);}
            else {m.component("comp1").physics("emw").create("impTi","Impedance",2);PhysicsFeature imp=m.component("comp1").physics("emw").feature("impTi");imp.selection().set(metal);imp.set("surfaceImpedanceDefinition","userdefined");imp.set("surfaceImpedance","Zs");}
            m.study().create("std1");m.study("std1").create("freq","Frequency");m.study("std1").feature("freq").set("plist","f0");

            progress("GENERATING_COMPLETE_PORT_BASIS");
            List<Map<String,Object>>expected=expectedOrders(lambda,requestedK[0],requestedK[1]);int expectedSpatial=expected.size(),expectedChannels=2*expectedSpatial;
            double[]paddingIncident=requestedIncidentK(lambda/1.10,theta,phi);List<Map<String,Object>>padding=expectedOrders(lambda/1.10,paddingIncident[0],paddingIncident[1]);
            Set<String>expectedPhysicalSet=independentOrderKeys(expected),expectedPaddingSet=independentOrderKeys(padding);
            int expectedPaddedChannels=authorizedPaddingChannels(id);if(!expectedPaddingSet.containsAll(expectedPhysicalSet)||2*padding.size()!=expectedPaddedChannels)throw new IllegalStateException("Independent retained padding contract changed");
            // A finite1.20f0 generator superset repairs COMSOL's observed label-range
            // truncation. Retain EXACTLY the original1.10f0 independently enumerated
            // channel set; neither physical study frequency nor padding count changes.
            port.set("MaximumFrequency","UserDefined");port.set("fmax","1.20*f0");port.runCommand("addDiffractionOrders");
            port=m.component("comp1").physics("emw").feature("port1");Map<String,Map<String,String>>superset=orderDescriptors(port);List<Map<String,Object>>removedSuperset=new ArrayList<>();
            for(Map.Entry<String,Map<String,String>>e:superset.entrySet())if(!expectedPaddingSet.contains(nativeDescriptorCanonicalKey(e.getValue()))){Map<String,String>d=e.getValue();removedSuperset.add(obj("key",e.getKey(),"feature_tag",d.get("feature_tag"),"port_name",d.get("port_name"),"native_m",Integer.parseInt(d.get("m")),"native_n",Integer.parseInt(d.get("n"))));m.component("comp1").physics("emw").feature("port1").feature().remove(d.get("feature_tag"));}
            port=m.component("comp1").physics("emw").feature("port1");port.set("fmax","f0");Map<String,Map<String,String>>channels=orderDescriptors(port);
            if(!expectedPairsExactly(channels,expectedPaddingSet))throw new IllegalStateException("Single1.20f0 superset did not retain every exact1.10f0 independent padding channel pair; no solve");
            Map<String,Map<String,String>>physicalMutable=new LinkedHashMap<>();for(Map.Entry<String,Map<String,String>>e:channels.entrySet())if(expectedPhysicalSet.contains(nativeDescriptorCanonicalKey(e.getValue())))physicalMutable.put(e.getKey(),e.getValue());
            Map<String,Map<String,String>>physical=Collections.unmodifiableMap(physicalMutable);Set<String>physicalKeys=Collections.unmodifiableSet(new LinkedHashSet<>(physical.keySet()));
            if(!expectedPairsExactly(physical,expectedPhysicalSet)||physical.size()+1!=expectedChannels)throw new IllegalStateException("Physical basis is not exact independent shifted lattice x2");
            List<Map<String,Object>>basis=new ArrayList<>();basis.add(obj("port_name","1","key","0,0,co","propagating",true));
            for(Map.Entry<String,Map<String,String>>e:channels.entrySet())basis.add(obj("port_name",e.getValue().get("port_name"),"key",e.getKey(),"propagating",physicalKeys.contains(e.getKey()),"feature_type",e.getValue().get("type")));
            emit("CODEX_PORT_BASIS_JSON",obj("case_id",caseId,"proof_schema",PROOF_SCHEMA,"wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"propagating_spatial_orders",expectedSpatial,"polarization_channels",expectedChannels,"total_port_channels_including_evanescent",channels.size()+1,"evanescent_channels",channels.size()-physical.size(),"padding_generation_frequency_ratio",1.10,"superset_generator_frequency_ratio",1.20,"independent_padding_set_exact",true,"single_finite_generation_attempt",true,"removed_superset_channels",removedSuperset,"physical_frequency_restored",true,"channels",basis,"independent_expected_lattice_orders",expected,"independent_retained_padding_orders",padding));
            progress("MESHING");m.component("comp1").mesh().create("mesh1");MeshSequence ms=m.component("comp1").mesh("mesh1");
            boolean short9=mesh.equals("m4_air28_film02_short9");double hb=2.8*UM,hs=(short9?.75:.875)*UM,he=(short9?.225:.2625)*UM;boolean localTop=false;
            // One bounded top-face refinement only; same geometry, film, Ti and edges.
            ms.feature("size").set("custom","on");ms.feature("size").set("hmax",hb);ms.feature("size").set("hmin",he/3);
            if(localTop){ms.create("sizeTop","Size");ms.feature("sizeTop").selection().geom("geom1",2);ms.feature("sizeTop").selection().set(top);ms.feature("sizeTop").set("custom","on");ms.feature("sizeTop").set("hmax",2.4*UM);ms.feature("sizeTop").set("hmin",he/3);}
            ms.create("sizeTi","Size");ms.feature("sizeTi").selection().geom("geom1",2);ms.feature("sizeTi").selection().set(metal);ms.feature("sizeTi").set("custom","on");ms.feature("sizeTi").set("hmax",hs);ms.feature("sizeTi").set("hmin",he/3);
            TreeSet<Integer> metalEdges=new TreeSet<>();for(int f:metal)for(int e:gi.getAdj(2,1,f))metalEdges.add(e);
            ms.create("sizeEdge","Size");ms.feature("sizeEdge").selection().geom("geom1",1);ms.feature("sizeEdge").selection().set(ints(metalEdges));ms.feature("sizeEdge").set("custom","on");ms.feature("sizeEdge").set("hmax",he);ms.feature("sizeEdge").set("hmin",he/3);
            double hFilm=(short9?.2:.5)*UM;
            ms.create("sizeFilm","Size");ms.feature("sizeFilm").selection().geom("geom1",3);ms.feature("sizeFilm").selection().set(new int[]{filmDomain});ms.feature("sizeFilm").set("custom","on");ms.feature("sizeFilm").set("hmax",hFilm);ms.feature("sizeFilm").set("hmin",he/3);
            for(int j=0;j<3;j++){String t="id"+j;ms.create(t,"IdenticalMesh");ms.feature(t).selection("group1").set(new int[]{sides[j]});ms.feature(t).selection("group2").set(new int[]{sides[j+3]});}
            ms.create("ftet1","FreeTet");ms.run();
            double[][] xyz=ms.getVertex();int[][] tri=ms.getElem("tri");int[] triEntities=ms.getElemEntity("tri");boolean meshMatch=true;List<Map<String,Object>> matches=new ArrayList<>();
            for(int j=0;j<3;j++){double ang=j*Math.PI/3,dx=-PITCH*Math.cos(ang),dy=-PITCH*Math.sin(ang);SortedMap<String,Integer>a=faceTriangles(xyz,tri,triEntities,sides[j],dx,dy),b=faceTriangles(xyz,tri,triEntities,sides[j+3],0,0);boolean eq=!a.isEmpty()&&a.equals(b);meshMatch&=eq;matches.add(obj("source_face",sides[j],"destination_face",sides[j+3],"translation_m",new double[]{dx,dy,0},"source_triangles",a.size(),"destination_triangles",b.size(),"canonical_triangle_multiset_equal",eq));}
            if(!meshMatch)throw new IllegalStateException("Actual paired triangle meshes not translation identical");
            progress("COMPILING_EQUATIONS_ONLY");m.study("std1").createAutoSequences("all");String solver=null,varTag=null;List<Map<String,Object>> solverAudit=new ArrayList<>();
            for(String s:m.sol().tags()){String found=null;for(String t:m.sol(s).feature().tags()){String type=m.sol(s).feature(t).getType();solverAudit.add(obj("solver",s,"feature",t,"type",type));if(type.equals("Variables"))found=t;}if(found!=null){if(solver!=null)throw new IllegalStateException("Ambiguous automatic solver sequences");solver=s;varTag=found;}}
            if(solver==null)throw new IllegalStateException("No actual Variables feature for presolve DOF inventory");
            m.sol(solver).run(varTag);int dof=m.sol(solver).feature(varTag).xmeshInfo().nDofs();int tets=ms.getNumElem("tet");double minQuality=ms.getMinQuality("tet");
            double[] presolveIncidentK={global(m,"real(emw.kIncx_1)"),global(m,"real(emw.kIncy_1)"),global(m,"real(emw.kIncz_1)")};
            boolean requestedIncidence=true;for(int j=0;j<3;j++)requestedIncidence&=Double.isFinite(presolveIncidentK[j])&&Math.abs(presolveIncidentK[j]-requestedK[j])<.1;
            boolean floquetBinding=true;List<Map<String,Object>> floquetReadback=new ArrayList<>();
            for(int j=0;j<3;j++){PhysicsFeature pc=m.component("comp1").physics("emw").feature("pc"+j);boolean ok="Floquet".equals(pc.getString("PeriodicType"))&&"FromPeriodicPort".equals(pc.getString("Floquet_source"));floquetBinding&=ok;floquetReadback.add(obj("tag","pc"+j,"PeriodicType",pc.getString("PeriodicType"),"Floquet_source",pc.getString("Floquet_source"),"verified",ok));}
            if(!requestedIncidence||!floquetBinding)throw new IllegalStateException("Requested global incidence or Floquet source was not reproduced by compiled equations");
            Map<String,Object>presolvePortProof=presolvePortWavevectorProof(m,top,topArea,lambda,theta,phi,requestedK,presolveIncidentK,expected,padding,channels,physicalKeys,removedSuperset);
            Map<String,Object> inventory=obj("case_id",caseId,"proof_schema",PROOF_SCHEMA,"wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"polarization",pol,"expected_requested_k_in_rad_per_m",requestedK,"actual_k_in_rad_per_m",presolveIncidentK,"requested_incidence_verified",requestedIncidence,"global_phi_verified",requestedIncidence,"side_floquet_source_verified",floquetBinding,"side_floquet_readback",floquetReadback,"state","MESH_PREFLIGHT_ONLY","solver_executed",false,"geometry_verified",geometryOK,"periodic_mesh_verified",meshMatch,"mesh",mesh,"propagating_spatial_orders",expectedSpatial,"polarization_channels",expectedChannels,"total_port_channels_including_evanescent",channels.size()+1,"tetrahedra",tets,"degrees_of_freedom",dof,"reported_dof_kind","actual","min_tetrahedron_quality",minQuality,"mesh_volume_m3",ms.getVolume("tet"),"geometry_volume_m3",volume,"periodic_mesh_pairs",matches,"solver_features",solverAudit,"second_order_curl_elements_requested",true,"finished_utc",Instant.now().toString());
            inventory.put("electric_curl_element_order",fieldOrder);inventory.put("height_um",HEIGHT/UM);inventory.put("air_cap_um",cap/UM);inventory.put("top_plane_z_um",ztop/UM);inventory.put("fixed_uniform_air_source_plane_z_um",sourcePlaneZUm);inventory.put("presolve_port_wavevector_proof",presolvePortProof);
            String meshSha=meshFingerprint(ms);inventory.put("mesh_fingerprint_sha256",meshSha);inventory.put("mesh_fingerprint_scheme","SHA256_RAW_COORDINATE_BITS_TET_TRI_CONNECTIVITY_ENTITY_IDS_V1");inventory.put("mesh_hAir_um",hb/UM);inventory.put("mesh_top_port_hmax_um",localTop?2.4:null);inventory.put("mesh_hLocalTop_um",localTop?2.4:null);inventory.put("mesh_local_variant",short9?"short9_film02_Ti075_edge0225":"original_m4_air28");
            boolean topSizeVerified=localTop?Math.abs(ms.feature("sizeTop").getDouble("hmax")-2.4*UM)<1e-18:!Arrays.asList(ms.feature().tags()).contains("sizeTop");inventory.put("mesh_local_top_boundary_size_readback_verified",topSizeVerified);if(!topSizeVerified)throw new IllegalStateException("Id-bound localtopmeshsize readback mismatch");inventory.put("mesh_hTi_um",hs/UM);inventory.put("mesh_hEdge_um",he/UM);inventory.put("mesh_hFilm_um",hFilm/UM);inventory.put("film_hmax_um",hFilm/UM);inventory.put("film_volume_m3",filmVolume);inventory.put("film_domain",filmDomain);inventory.put("air_domain",airDomain);inventory.put("floor_fill_thickness_um",filmUm);inventory.put("minimum_quality_admission",.20);
            double actualAir=ms.feature("size").getDouble("hmax"),actualTi=ms.feature("sizeTi").getDouble("hmax"),actualEdge=ms.feature("sizeEdge").getDouble("hmax"),actualFilm=ms.feature("sizeFilm").getDouble("hmax");
            boolean exactSizeReadback=Math.abs(actualAir-hb)<1e-18&&Math.abs(actualTi-hs)<1e-18&&Math.abs(actualEdge-he)<1e-18&&Math.abs(actualFilm-hFilm)<1e-18;
            double filmIntensityDecayUm=lambda/(4*Math.PI*kfilm),tiFieldSkinUm=lambda/(2*Math.PI*kti);boolean lossMeshAdmitted=!short9||(lambda==9&&actualFilm/UM<=.5*filmIntensityDecayUm);
            if(!exactSizeReadback||!lossMeshAdmitted)throw new IllegalStateException("Exact case-bound mesh/readback or conditional9 attenuation-length admission failed");
            Map<String,Object>lossMeshProof=obj("profile",mesh,"conditional9",short9,"actual_profile_size_readback_verified",exactSizeReadback,"actual_air_hmax_um",actualAir/UM,"actual_Ti_boundary_hmax_um",actualTi/UM,"actual_edge_hmax_um",actualEdge/UM,"actual_film_hmax_um",actualFilm/UM,"film_intensity_decay_length_um",filmIntensityDecayUm,"film_hmax_over_intensity_decay",actualFilm/UM/filmIntensityDecayUm,"conditional9_maximum_film_hmax_over_intensity_decay",.5,"Ti_field_skin_depth_um",tiFieldSkinUm,"Ti_model","local_impedance_boundary_no_Ti_volume_skin_mesh","conditional_loss_mesh_admission",lossMeshAdmitted,"mesh_convergence_established",false);
            inventory.put("loss_region_mesh_admission",lossMeshProof);inventory.put("case_bound_mesh_size_readback_verified",exactSizeReadback);
            int scopedDofLimit=8000000; // Exactly eleven id-bound cases; every case preserves the original coarse mesh.
            boolean feasible=dof>0&&dof<=scopedDofLimit&&tets>0&&tets<=2500000&&Double.isFinite(minQuality)&&minQuality>=.20;inventory.put("eligible_for_manual_review",feasible);
            Map<String,Object> linearOverride=exactOutOfCore(m,solver,true);emit("CODEX_LINEAR_SOLVER_OVERRIDE_JSON",linearOverride);inventory.put("linear_solver_override",linearOverride);
            m.save(root+"\\mesh_only.mph");emit("CODEX_MESH_JSON",inventory);
            if(meshOnly){emit("CODEX_DONE_JSON",obj("mesh_only",true,"solver_executed",false,"model_file",root+"\\mesh_only.mph","eligible_for_manual_review",feasible));progress("MESH_PREFLIGHT_ONLY");return;}
            if(!feasible||reviewedDof!=dof||reviewedTet!=tets||reviewedMeshSha==null||!reviewedMeshSha.matches("[a-f0-9]{64}")||!reviewedMeshSha.equals(meshSha))throw new IllegalStateException("Resource/quality or exact reviewed mesh SHA/DOF/tet admission mismatch; no physical solve");
            exactOutOfCore(m,solver,false);
            progress("SOLVING");long begin=System.nanoTime();m.sol(solver).runAll();solverCompleted=true;double seconds=(System.nanoTime()-begin)*1e-9;
            progress("SAVING_SOLVED_UNEVALUATED");m.save(root+"\\solved_unevaluated.mph");solvedArtifactSaved=true;
            emit("CODEX_SOLVED_ARTIFACT_JSON",obj("case_id",caseId,"state","SOLVED_UNEVALUATED_SAVED","solver_executed",true,
                "qualified",false,"model_path",root+"\\solved_unevaluated.mph","solve_seconds",seconds,"saved_utc",Instant.now().toString()));
            progress("EVALUATING");evaluationCheckpoint("EVALUATING_TOP_NET_AND_SOURCE_PORT",obj());
            String pout="(emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz)";
            double topnet=integral(m,top,pout),Aflux=-topnet,Rnet=1+topnet,Rs=global(m,"abs(emw.S11)^2"),Rtotal=Rs,nonzero=0;
            evaluationCheckpoint("TOP_NET_AND_SOURCE_PORT_EVALUATED",obj("A_flux",Aflux,"R_total_net_flux",Rnet,"R_zeroth_co",Rs));
            double k0=2*Math.PI/(lambda*UM);int[] actualSpatialCounts=new int[expected.size()];int incidentExpectedIndex=expectedIndex(expected,requestedK[0],requestedK[1]);if(incidentExpectedIndex<0)throw new IllegalStateException("Missing requested zero order in shifted lattice");actualSpatialCounts[incidentExpectedIndex]=1;
            boolean actualWavesOK=true;
            List<Map<String,Object>> powers=new ArrayList<>();powers.add(obj("port_name","1","key","0,0,co","reflected_power_over_Pin",Rs));
            boolean orderFinite=Double.isFinite(Rs);for(Map.Entry<String,Map<String,String>> e:channels.entrySet()){
                String pn=e.getValue().get("port_name");boolean prop=physicalKeys.contains(e.getKey());
                evaluationCheckpoint("EVALUATING_ORDER",obj("port_name",pn,"order_key",e.getKey(),"propagating",prop));
                double kx=avg(m,top,topArea,"real(emw.kModex_"+pn+")"),ky=avg(m,top,topArea,"real(emw.kModey_"+pn+")"),kz=avg(m,top,topArea,"real(emw.kModez_"+pn+")"),ikz=avg(m,top,topArea,"imag(emw.kModez_"+pn+")");
                int expectedId=expectedIndex(expected,kx,ky);boolean actualProp=kx*kx+ky*ky<k0*k0*(1-1e-10);
                actualWavesOK&=Double.isFinite(kx)&&Double.isFinite(ky)&&Double.isFinite(kz)&&Double.isFinite(ikz)&&actualProp==prop;
                if(prop){actualWavesOK&=expectedId>=0&&Math.abs(ikz)<.1&&Math.abs((kx*kx+ky*ky+kz*kz)/(k0*k0)-1)<1e-7;if(expectedId>=0)actualSpatialCounts[expectedId]++;}
                Map<String,Object> row=obj("port_name",pn,"key",e.getKey(),"propagating",prop,"actual_wavevector_propagating",actualProp,"independent_expected_spatial_index",expectedId,"kx_rad_per_m",kx,"ky_rad_per_m",ky,"real_kz_rad_per_m",kz,"imag_kz_rad_per_m",ikz);
                if(prop){String svar=scatteringVariable(m,pn);double p=global(m,"abs("+svar+")^2");Rtotal+=p;if(!e.getKey().equals("0,0,orthogonal"))nonzero+=p;row.put("S_variable",svar);row.put("reflected_power_over_Pin",p);orderFinite&=Double.isFinite(p);}powers.add(row);
                evaluationCheckpoint("ORDER_EVALUATED",row);
            }
            for(int n:actualSpatialCounts)actualWavesOK&=n==2;
            evaluationCheckpoint("ORDERS_EVALUATED",obj("R_total",Rtotal,"R_nonzero_orders",nonzero,"actual_wavevectors_match_independent_lattice",actualWavesOK));
            evaluationCheckpoint("EVALUATING_INCIDENT_WAVE_AND_ZERO_ORDER_PROJECTION",obj());
            double[] incidentK=new double[]{avg(m,top,topArea,"real(emw.kIncx_1)"),avg(m,top,topArea,"real(emw.kIncy_1)"),avg(m,top,topArea,"real(emw.kIncz_1)")};
            boolean incidentDirection=true;for(int j=0;j<3;j++)incidentDirection&=Double.isFinite(incidentK[j])&&Math.abs(incidentK[j]-requestedK[j])<.1;
            String topDephase="exp(i*(("+incidentK[0]+"[1/m])*x+("+incidentK[1]+"[1/m])*y))";
            double Z0=global(m,"Z0_const"),ct=Math.cos(tr);double erx=avg(m,top,topArea,"real(emw.Ex*("+topDephase+"))"),eix=avg(m,top,topArea,"imag(emw.Ex*("+topDephase+"))"),ery=avg(m,top,topArea,"real(emw.Ey*("+topDephase+"))"),eiy=avg(m,top,topArea,"imag(emw.Ey*("+topDephase+"))"),hrx=avg(m,top,topArea,"real(emw.Hx*("+topDephase+"))"),hix=avg(m,top,topArea,"imag(emw.Hx*("+topDephase+"))"),hry=avg(m,top,topArea,"real(emw.Hy*("+topDephase+"))"),hiy=avg(m,top,topArea,"imag(emw.Hy*("+topDephase+"))");
            double cp=Math.cos(Math.toRadians(phi)),sp=Math.sin(Math.toRadians(phi));
            double epr=cp*erx+sp*ery,epi=cp*eix+sp*eiy,esr=-sp*erx+cp*ery,esi=-sp*eix+cp*eiy;
            double hpr=cp*hrx+sp*hry,hpi=cp*hix+sp*hiy,hsr=-sp*hrx+cp*hry,hsi=-sp*hix+cp*hiy;
            double PinField=topArea*(sq(epr-Z0*ct*hsr,epi-Z0*ct*hsi)/ct+ct*sq(esr+Z0*hpr/ct,esi+Z0*hpi/ct))/(8*Z0);
            double pref0=topArea*(sq(epr+Z0*ct*hsr,epi+Z0*ct*hsi)/ct+ct*sq(esr-Z0*hpr/ct,esi-Z0*hpi/ct))/(8*Z0);
            evaluationCheckpoint("INCIDENT_WAVE_AND_ZERO_ORDER_PROJECTION_EVALUATED",obj("incident_k_rad_per_m",incidentK,"incident_wavevector_points_down",incidentDirection,"Pin_top_EH_zeroth_projection_W",PinField,"Pref_top_EH_zeroth_projection_W",pref0,"role","Original top raw-H diagnostic, not primary incident-power admission"));
            evaluationCheckpoint("EVALUATING_NORMALIZED_SOURCE_MODE_DIRECTION",obj());
            double[] modeComponents={integral(m,top,"abs(emw.tEmodex_1)^2"),integral(m,top,"abs(emw.tEmodey_1)^2"),integral(m,top,"abs(emw.tEmodez_1)^2")};
            double modeAxisX=pol.equals("s")?-sp:cp,modeAxisY=pol.equals("s")?cp:sp;
            double modeTotal=modeComponents[0]+modeComponents[1]+modeComponents[2],targetModeFraction=integral(m,top,"abs(("+modeAxisX+")*emw.tEmodex_1+("+modeAxisY+")*emw.tEmodey_1)^2")/modeTotal,modeLeakage=1-targetModeFraction;
            double[] modeFractions={modeComponents[0]/modeTotal,modeComponents[1]/modeTotal,modeComponents[2]/modeTotal};
            boolean modeDirection=Double.isFinite(modeTotal)&&modeTotal>0&&Double.isFinite(modeLeakage)&&targetModeFraction>=1-MODE_DIRECTION_LEAKAGE_MAX&&modeLeakage<=MODE_DIRECTION_LEAKAGE_MAX;
            evaluationCheckpoint("NORMALIZED_SOURCE_MODE_DIRECTION_EVALUATED",obj("normalized_tangential_mode_E_component_fractions_xyz",modeFractions,"target_global_component_fraction",targetModeFraction,"orthogonal_leakage_fraction",modeLeakage,"mode_direction_passed",modeDirection));
            Map<String,Object> powerFit=internalElectricPowerFit(m,topArea,zland,cap,lambda,incidentK);
            double PinIndependent=((Number)powerFit.get("Pin_W")).doubleValue();
            evaluationCheckpoint("EVALUATING_SURFACE_DISSIPATION",obj());
            String q="0.5*real(Zs)*(abs(emw.Hx)^2+abs(emw.Hy)^2+abs(emw.Hz)^2-abs(nx*emw.Hx+ny*emw.Hy+nz*emw.Hz)^2)";
            double afloor=pec?0:integral(m,floor,q),awalls=pec?0:integral(m,walls,q),aland=pec?0:integral(m,land,q),Asurf=afloor+awalls+aland,edgeLoss=pec?0:integral(m,metal,"("+q+")*("+edgeMask(flat,HEIGHT)+")");
            evaluationCheckpoint("SURFACE_DISSIPATION_EVALUATED",obj("A_floor",afloor,"A_walls",awalls,"A_land",aland,"A_surface_Ti_H",Asurf,"edge_associated_absorption",edgeLoss));
            evaluationCheckpoint("EVALUATING_EXPLICIT_FILM_AND_TI_E_DISSIPATION",obj());
            String e2="abs(emw.Ex)^2+abs(emw.Ey)^2+abs(emw.Ez)^2";
            double Afilm=volumeIntegral(m,new int[]{filmDomain},"-0.5*(2*pi*f0)*epsilon0_const*imag(eps_f)*("+e2+")"),AfilmQh=volumeIntegral(m,new int[]{filmDomain},"emw.Qh"),AtiE=integral(m,metal,"0.5*real(1/Zs)*("+e2+"-abs(nx*emw.Ex+ny*emw.Ey+nz*emw.Ez)^2)"),totalLoss=Asurf+Afilm;
            double filmTiH=integral(m,filmTi,q),airTiH=integral(m,airTi,q);
            evaluationCheckpoint("FILM_AND_TI_DISSIPATION_EVALUATED",obj("A_film",Afilm,"A_film_builtin_Qh",AfilmQh,"A_Ti_H",Asurf,"A_Ti_E",AtiE,"A_loss_total",totalLoss,"A_Ti_film_side",filmTiH,"A_Ti_air_side",airTiH));
            evaluationCheckpoint("EVALUATING_PERIODIC_SIDE_POWER",obj());
            double sideMismatch=0;List<Map<String,Object>> sidePower=new ArrayList<>();for(int j=0;j<3;j++){double p=integral(m,new int[]{sides[j]},pout),q2=integral(m,new int[]{sides[j+3]},pout);sideMismatch=Math.max(sideMismatch,Math.abs(p+q2));sidePower.add(obj("source_W",p,"destination_W",q2,"sum_W",p+q2));}
            double Ribc=global(m,"R_ibc"),aibc=1-Ribc,pinerr=Math.abs(PinIndependent-1),rawTopPinerr=Math.abs(PinField-1),orderdiff=Math.abs(Rtotal-Rnet),absclose=Math.abs(Aflux-totalLoss),relclose=absclose/Math.max(Math.abs(Aflux),1e-12);
            evaluationCheckpoint("ASSEMBLING_FINAL_CHECKS",obj("R_total",Rtotal,"A_flux",Aflux,"A_surface_Ti_H",Asurf,"incident_power_relative_error",pinerr,"R_orders_vs_net_flux_abs_difference",orderdiff,"flux_total_loss_absolute_mismatch",absclose,"flux_total_loss_relative_mismatch",relclose,"side_pair_power_mismatch_over_Pin",sideMismatch));
            double filmQhAbs=Math.abs(Afilm-AfilmQh),filmQhRel=filmQhAbs/Math.max(Math.abs(Afilm),1e-12),TiHEAbs=Math.abs(Asurf-AtiE),TiHERel=TiHEAbs/Math.max(Math.abs(Asurf),1e-12),portLossAbs=Math.abs(1-Rtotal-totalLoss),portLossRel=portLossAbs/Math.max(Math.abs(1-Rtotal),1e-12);
            Map<String,Object> strictDiagnostics=obj("absolute_threshold",STRICT_ABS_DIAGNOSTIC,"R_orders_vs_net_flux",orderdiff<=STRICT_ABS_DIAGNOSTIC,"flux_vs_total_loss_absolute",absclose<=STRICT_ABS_DIAGNOSTIC,"film_manual_vs_builtin_Qh_absolute",filmQhAbs<=STRICT_ABS_DIAGNOSTIC,"Ti_H_vs_E_loss_absolute",TiHEAbs<=STRICT_ABS_DIAGNOSTIC,"port_vs_total_loss_absolute",portLossAbs<=STRICT_ABS_DIAGNOSTIC,"role","Retained strict absolute0.002 diagnostics, not exploratory admission gates");
            Map<String,Object> checks=obj("geometry_verified",geometryOK,"periodic_mesh_verified",meshMatch,"all_propagating_channels_present",physicalKeys.size()+1==expectedChannels&&validPhysicalFeaturePairs(physical,expectedSpatial),"actual_wavevectors_match_independent_lattice",actualWavesOK,"incident_wavevector_points_down",incidentDirection,"finite",orderFinite&&Double.isFinite(Aflux)&&Double.isFinite(Asurf)&&Double.isFinite(AtiE)&&Double.isFinite(Afilm)&&Double.isFinite(AfilmQh)&&Double.isFinite(PinField),"incident_power",pinerr<=.005,"R_orders_vs_net_flux",orderdiff<=ABS_POWER_GATE,"periodic_side_power",sideMismatch<=.002,"flux_vs_total_loss_absolute",absclose<=ABS_POWER_GATE,"flux_vs_total_loss_relative",relclose<=REL_LOSS_GATE,"passive",Aflux>=-.002&&Aflux<=1.002&&Rtotal>=-.002&&Rtotal<=1.002);
            checks.put("requested_global_incidence",incidentDirection&&requestedIncidence);checks.put("side_floquet_source_binding",floquetBinding);
            checks.put("explicit_global_source_binding",sourceBinding);checks.put("normalized_source_mode_global_direction",modeDirection);
            checks.put("internal_E_fit_finite",Double.isFinite(PinIndependent));checks.put("internal_E_fit_condition",Boolean.TRUE.equals(powerFit.get("condition_passed")));
            checks.put("internal_E_fit_relative_L2_residual",Boolean.TRUE.equals(powerFit.get("residual_passed")));checks.put("internal_E_fit_all_LOO_Pin",Boolean.TRUE.equals(powerFit.get("LOO_Pin_passed")));checks.put("internal_E_fit_cut_areas",Boolean.TRUE.equals(powerFit.get("all_cut_areas_match_geometry")));
            checks.put("film_parameter_binding",filmParameterBinding);checks.put("explicit_film_passive",Double.isFinite(Afilm)&&Afilm>=-1e-10);
            checks.put("film_manual_vs_builtin_Qh_absolute",filmQhAbs<=ABS_POWER_GATE);checks.put("film_manual_vs_builtin_Qh_relative",filmQhRel<=REL_LOSS_GATE);
            checks.put("Ti_H_vs_E_loss_absolute",TiHEAbs<=ABS_POWER_GATE);checks.put("Ti_H_vs_E_loss_relative",TiHERel<=REL_LOSS_GATE);
            checks.put("port_vs_total_loss_absolute",portLossAbs<=ABS_POWER_GATE);checks.put("port_vs_total_loss_relative",portLossRel<=REL_LOSS_GATE);
            checks.put("loss_partition_complete",Math.abs(filmTiH+airTiH-Asurf)<=1e-8);checks.put("total_loss_passive",totalLoss>=-.002&&totalLoss<=1.002);
            boolean passed=true;for(Object v:checks.values())passed&=Boolean.TRUE.equals(v);String state=passed?"PASSED_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY":"FAILED_FLOOR_FILL_CAVITY_DIAGNOSTIC";
            Map<String,Object> metrics=obj("case_id",caseId,"state",state,"solver_executed",true,"geometry",geometry,"surface",surface,"wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"polarization",pol,"mesh",mesh,"solve_seconds",seconds,"R_total",Rtotal,"R_total_net_flux",Rnet,"R_zeroth_co",Rs,"R_nonzero_orders",nonzero,"R_zeroth_projection",pref0/PinField,"A_flux",Aflux,"A_surface_Ti_H",Asurf,"A_floor",afloor,"A_walls",awalls,"A_land",aland,"T",0,"Pin_specified_W",1,"Pin_field_zeroth_projection_W",PinField,"incident_power_relative_error",pinerr,"R_orders_vs_net_flux_abs_difference",orderdiff,"side_pair_power_mismatch_over_Pin",sideMismatch,"flux_total_loss_absolute_mismatch",absclose,"flux_total_loss_relative_mismatch",relclose,"R_analytic_ibc",Ribc,"numerical_A_relative_error",Math.abs(Aflux-aibc)/aibc,"Pout_top_net_W",topnet,"side_pair_powers",sidePower,"edge_strip_width_um",.5,"edge_associated_absorption",edgeLoss,"edge_associated_loss_fraction",Asurf>0?edgeLoss/Asurf:0,"edge_fraction_requires_material_model_review",Asurf>0&&edgeLoss/Asurf>.05,"incident_k_rad_per_m",new double[]{avg(m,top,topArea,"real(emw.kIncx_1)"),avg(m,top,topArea,"real(emw.kIncy_1)"),avg(m,top,topArea,"real(emw.kIncz_1)")},"top_zero_order_E_real",new double[]{erx,ery,0},"top_zero_order_E_imag",new double[]{eix,eiy,0},"top_zero_order_H_real",new double[]{hrx,hry,0},"top_zero_order_H_imag",new double[]{hix,hiy,0},"order_powers",powers,"checks",checks,"scope","A3 H="+HEIGHT/UM+"um lower0.."+filmUm+"um explicit SiO2 floor-fill only, bare upper walls/land,"+lambda+"um theta="+theta+"deg phi="+phi+"deg "+pol+", local Ti IBC; partial angular exploration NOT fullhemisphere/completeband/measuredTA2/formalqualification","started_utc",started,"finished_utc",Instant.now().toString());
            metrics.put("numerical_A_relative_error",null);metrics.put("height_um",HEIGHT/UM);metrics.put("air_cap_um",cap/UM);metrics.put("top_plane_z_um",ztop/UM);metrics.put("fixed_uniform_air_source_plane_z_um",sourcePlaneZUm);metrics.put("nfilm",nfilm);metrics.put("kfilm",kfilm);metrics.put("floor_fill_thickness_um",filmUm);metrics.put("film_volume_m3",filmVolume);metrics.put("film_hmax_um",hFilm/UM);metrics.put("not_conformal_wall_or_land_coating",true);metrics.put("A_film",Afilm);metrics.put("A_film_builtin_Qh",AfilmQh);metrics.put("A_Ti_H",Asurf);metrics.put("A_Ti_E",AtiE);metrics.put("A_loss_total",totalLoss);metrics.put("A_Ti_film_side",filmTiH);metrics.put("A_Ti_air_side",airTiH);metrics.put("mesh_fingerprint_sha256",meshSha);metrics.put("reviewed_mesh_fingerprint_sha256",reviewedMeshSha);metrics.put("film_manual_Qh_absolute_mismatch",Math.abs(Afilm-AfilmQh));metrics.put("film_manual_Qh_relative_mismatch",Math.abs(Afilm-AfilmQh)/Math.max(Math.abs(Afilm),1e-12));metrics.put("Ti_H_E_absolute_mismatch",Math.abs(Asurf-AtiE));metrics.put("Ti_H_E_relative_mismatch",Math.abs(Asurf-AtiE)/Math.max(Math.abs(Asurf),1e-12));metrics.put("planar_TMM_is_not_cavity_acceptance_reference",true);metrics.put("edge_loss_scope","Original18Ti sharp-edge segments only; extra film-termination contact-line loss is not isolated");
            metrics.put("A_port",1-Rtotal);metrics.put("A_loss_total_Ti_E",Afilm+AtiE);metrics.put("port_loss_absolute_mismatch",portLossAbs);metrics.put("port_loss_relative_mismatch",portLossRel);metrics.put("strict_0p002_diagnostics",strictDiagnostics);
            metrics.put("solver_completed",true);metrics.put("nti",nti);metrics.put("kti",kti);metrics.put("titanium_parameter_binding_verified",titaniumParameterBinding);metrics.put("mesh_hLocalTop_um",localTop?2.4:null);metrics.put("mesh_local_top_boundary_size_readback_verified",topSizeVerified);metrics.put("eps_f_real",nfilm*nfilm-kfilm*kfilm);metrics.put("eps_f_imag",-2*nfilm*kfilm);metrics.put("film_parameter_binding_verified",filmParameterBinding);
            metrics.put("proof_schema",PROOF_SCHEMA);metrics.put("requested_k_in_rad_per_m",requestedK);metrics.put("incoming_E_unit_global",incomingE);metrics.put("outgoing_E_unit_global",outgoingE);metrics.put("port_azimuth_offset_deg",PORT_AZIMUTH_OFFSET_DEG);metrics.put("requested_incidence_verified",incidentDirection&&requestedIncidence);metrics.put("global_phi_verified",incidentDirection&&requestedIncidence);metrics.put("side_floquet_source_verified",floquetBinding);metrics.put("side_floquet_readback",floquetReadback);metrics.put("propagating_spatial_orders",expectedSpatial);metrics.put("polarization_channels",expectedChannels);metrics.put("total_port_channels_including_evanescent",channels.size()+1);
            metrics.put("linear_solver_override",linearOverride);metrics.put("linear_solver_override_reverified_before_physics",true);metrics.put("normal_incidence_polarization_definition",theta==0?(pol.equals("s")?"global_Ey":"global_Ex"):"not_normal_incidence");metrics.put("polarization_definition",pol.equals("s")?"true global s: incoming/outgoing=(-sin(phi),cos(phi),0); tangential mode projected onto(-sin(phi),cos(phi),0)":"true global p: incoming=(cos(theta)cos(phi),cos(theta)sin(phi),+sin(theta)); outgoing=(cos(theta)cos(phi),cos(theta)sin(phi),-sin(theta)); tangential mode projected onto(cos(phi),sin(phi),0)");
            metrics.put("scope","A3 H="+HEIGHT/UM+"um lower0.."+filmUm+"um explicit SiO2 floor-fill only, bare upper walls/land,"+lambda+"um theta="+theta+"deg global phi="+phi+"deg "+pol+", local Ti IBC; single-point exploration NOT conformal coating, measured TA2, full broadband/angular emissivity, mesh-convergence or edge qualification");
            metrics.put("presolve_port_wavevector_proof",presolvePortProof);
            metrics.put("degrees_of_freedom",dof);metrics.put("tetrahedra",tets);metrics.put("min_tetrahedron_quality",minQuality);metrics.put("mesh_hAir_um",hb/UM);metrics.put("mesh_hTi_um",hs/UM);metrics.put("mesh_hEdge_um",he/UM);metrics.put("mesh_hFilm_um",hFilm/UM);
            metrics.put("exploratory_acceptance",obj("absolute_power_difference",ABS_POWER_GATE,"relative_absorption_difference",REL_LOSS_GATE,"periodic_side_power_absolute",.002,"passivity_absolute",.002,"relative_net_loss_denominator","abs(A_flux)","relative_port_loss_denominator","abs(1-R_total)","relative_Ti_H_E_denominator","abs(A_Ti_H)","relative_film_Qh_denominator","abs(A_film_manual)","strict_0p002_diagnostics_are_not_admission_gates",true));
            metrics.put("incident_k_rad_per_m",incidentK);
            metrics.put("actual_independent_spatial_channel_counts",actualSpatialCounts);
            metrics.put("electric_curl_element_order",fieldOrder);metrics.put("loss_region_mesh_admission",lossMeshProof);metrics.put("case_bound_mesh_size_readback_verified",exactSizeReadback);
            metrics.put("incident_power_method",INCIDENT_POWER_METHOD);metrics.put("Pin_independent_internal_E_fit_W",PinIndependent);metrics.put("internal_E_power_fit",powerFit);
            metrics.put("Pin_field_zeroth_projection_W",PinIndependent);metrics.put("Pref_field_zeroth_projection_W",powerFit.get("Pref_W"));
            metrics.put("R_zeroth_projection",((Number)powerFit.get("Pref_W")).doubleValue()/PinIndependent);
            metrics.put("Pin_top_EH_zeroth_projection_W",PinField);metrics.put("Pref_top_EH_zeroth_projection_W",pref0);
            Map<String,Object> metricSourceProof=new LinkedHashMap<>(sourceProof);metricSourceProof.put("requested",Arrays.asList(desiredGlobalE));metricSourceProof.put("verified",sourceBinding);metrics.put("explicit_global_source_binding",metricSourceProof);
            metrics.put("normalized_source_mode_E_component_fractions_xyz",modeFractions);metrics.put("normalized_source_mode_target_component_fraction",targetModeFraction);metrics.put("normalized_source_mode_target_tangential_unit_axis",new double[]{modeAxisX,modeAxisY,0});metrics.put("source_mode_projection_is_global_phi_aware",true);metrics.put("normalized_source_mode_orthogonal_leakage_fraction",modeLeakage);
            metrics.put("normalized_source_mode_component_scope","tangential_port_mode");
            metrics.put("original_top_raw_H_power_diagnostic",obj("Pin_W",PinField,"Pref_W",pref0,"R_zeroth_projection",pref0/PinField,"Pin_relative_error_vs_specified_1W",rawTopPinerr,"original_incident_power_threshold",.005,"original_incident_power_flag",rawTopPinerr<=.005,
                "role","Retained boundary derivative-field diagnostic; not the primary source-power admission measurement","saved_field_evidence","Original top curl-H trace differed from normalized-mode Maxwell H while source modal power and all-five-plane electric fit were consistent with specified 1W"));
            metrics.put("R_A_power_normalization","specified_1W_unchanged_not_fitted_Pin");metrics.put("measurement_definition_revision","Primary independent incident power uses five predetermined internal-cap Bloch-dephased E planes at fixed k0/kz with correct tangential p/s admittance; angle-correct raw-H top Pin/Pref and flag retained, not silently reclassified as passed");
            evaluationCheckpoint("SAVING_EVALUATED_MODEL",obj("numerical_checks_passed",passed));
            m.save(root+"\\model.mph");emit("CODEX_METRICS_JSON",metrics);emit("CODEX_DONE_JSON",obj("mesh_only",false,"solver_executed",true,"model_file",root+"\\model.mph","state",state,"numerical_checks_passed",passed));progress(state);
        } catch(Exception e){
            emit("CODEX_FAILURE_JSON",obj("case_id",caseId,"state","FAILED_EXCEPTION","failed_phase",lastPhase,
                "solver_completed",solverCompleted,"solved_artifact_saved",solvedArtifactSaved,
                "last_numerical_tag",lastNumericalTag,"last_numerical_kind",lastNumericalKind,"last_expression",lastExpression,
                "exception_type",e.getClass().getName(),"exception_message",e.getMessage(),"updated_utc",Instant.now().toString()));
            progress("FAILED_EXCEPTION");e.printStackTrace();throw e;
        }
        finally {if(m!=null)ModelUtil.remove("CavityFloorFilm");}
    }


    private static int ev=0;
    private static Map<String,Object> flatObj(Object... kv){Map<String,Object> r=new LinkedHashMap<>();for(int i=0;i<kv.length;i+=2)r.put((String)kv[i],kv[i+1]);return r;}
    private static String flatJson(Object x){
        if(x==null)return "null";if(x instanceof String)return "\""+((String)x).replace("\\","\\\\").replace("\"","\\\"").replace("\n","\\n").replace("\r","\\r")+"\"";
        if(x instanceof Number){if(!Double.isFinite(((Number)x).doubleValue()))throw new IllegalStateException("Nonfinite calibration metric");return x.toString();}
        if(x instanceof Boolean)return x.toString();StringJoiner j;
        if(x instanceof Map){j=new StringJoiner(",","{","}");for(Object k:((Map<?,?>)x).keySet())j.add(flatJson(k.toString())+":"+flatJson(((Map<?,?>)x).get(k)));return j.toString();}
        if(x instanceof Collection){j=new StringJoiner(",","[","]");for(Object q:(Collection<?>)x)j.add(flatJson(q));return j.toString();}
        if(x instanceof double[]){List<Double> r=new ArrayList<>();for(double q:(double[])x)r.add(q);return flatJson(r);}throw new IllegalArgumentException("Unsupported calibration JSON");
    }
    private static double integ(Model m,int face,String e){String tag="flatEv"+(++ev);m.result().numerical().create(tag,"IntSurface");try{NumericalFeature n=m.result().numerical(tag);n.selection().set(new int[]{face});n.set("data","dset1");n.set("expr",e);n.set("intorderactive","on");n.set("intorder",8);return n.getReal()[0][0];}finally{m.result().numerical().remove(tag);}}
    private static double glob(Model m,String e){String tag="flatEv"+(++ev);m.result().numerical().create(tag,"EvalGlobal");try{NumericalFeature n=m.result().numerical(tag);n.set("data","dset1");n.set("expr",e);return n.getReal()[0][0];}finally{m.result().numerical().remove(tag);}}
    public static void runPoint(String root,String fixedId,double lambdaIn,double thetaIn,double ntiIn,double ktiIn,String pol,double phi) throws Exception {
        System.out.println("CODEX_ARITHMETIC_SELF_TEST_JSON="+flatJson(FlatPrecheckEntry.pureArithmeticSelfTest()));System.out.flush();
        if(!"p".equals(pol)&&!"s".equals(pol))throw new IllegalArgumentException("Only explicit s/p flat source");
        boolean approved=false;
        if(fixedId.equals("F_ML6_A3H100_t2_l10p5_th0_phi0_s")&&lambdaIn==10.5&&thetaIn==0&&phi==0&&pol.equals("s")) approved=Math.abs(ntiIn-4.308219763636363)<1e-14&&Math.abs(ktiIn-20.63042690909091)<1e-14;
        if(fixedId.equals("F_ML7_A3H100_t2_l10p5_th0_phi0_p")&&lambdaIn==10.5&&thetaIn==0&&phi==0&&pol.equals("p")) approved=Math.abs(ntiIn-4.308219763636363)<1e-14&&Math.abs(ktiIn-20.63042690909091)<1e-14;
        if(fixedId.equals("F_ML8_A3H100_t2_l10p5_th30_phi30_s")&&lambdaIn==10.5&&thetaIn==30&&phi==30&&pol.equals("s")) approved=Math.abs(ntiIn-4.308219763636363)<1e-14&&Math.abs(ktiIn-20.63042690909091)<1e-14;
        if(fixedId.equals("F_ML9_A3H100_t2_l10p5_th30_phi30_p")&&lambdaIn==10.5&&thetaIn==30&&phi==30&&pol.equals("p")) approved=Math.abs(ntiIn-4.308219763636363)<1e-14&&Math.abs(ktiIn-20.63042690909091)<1e-14;
        if(!approved)throw new IllegalArgumentException("Flat ID/actual source/Ti outside closed depth-screen contract");
        double[][] points={{lambdaIn,thetaIn,ntiIn,ktiIn}};
        for(double[] point:points){
            double lambda=point[0],theta=point[1],tr=Math.toRadians(theta),ap=1e-6,H=10.5e-6,area=2*Math.sqrt(3)*ap*ap;
            String id=fixedId;Model m=null;ev=0;FlatPrecheckEntry.setProbeContext(id);
            try{
                m=ModelUtil.create("FlatAdmission");m.param().set("lambda0",lambda+"[um]");m.param().set("f0","c_const/lambda0");
                m.param().set("Zs","Z0_const/("+point[2]+"-i*"+point[3]+")");m.param().set("R_ibc",pol.equals("s")?"abs((Zs*cos("+theta+"[deg])-Z0_const)/(Zs*cos("+theta+"[deg])+Z0_const))^2":"abs((Zs-Z0_const*cos("+theta+"[deg]))/(Zs+Z0_const*cos("+theta+"[deg])))^2");
                m.component().create("comp1",true);m.component("comp1").geom().create("geom1",3);GeomSequence g=m.component("comp1").geom("geom1");g.lengthUnit("m");
                double[] x=new double[6],y=new double[6];for(int j=0;j<6;j++){double a=Math.PI/6+j*Math.PI/3;x[j]=ap/Math.cos(Math.PI/6)*Math.cos(a);y[j]=ap/Math.cos(Math.PI/6)*Math.sin(a);}
                g.create("wp","WorkPlane");g.feature("wp").geom().create("hex","Polygon");g.feature("wp").geom().feature("hex").set("x",x);g.feature("wp").geom().feature("hex").set("y",y);g.run("wp");g.create("ext","Extrude");g.feature("ext").set("workplane","wp");g.feature("ext").set("distance",H);g.run();
                double[][] v=g.getVertexCoord();int top=-1,bottom=-1,ref=-1;int[] sides=new int[6];
                for(int f=1;f<=g.getNFaces();f++){int[] ids=g.getAdj(2,0,f);boolean t=true,b=true;for(int z:ids){t&=Math.abs(v[2][z-1]-H)<1e-11;b&=Math.abs(v[2][z-1])<1e-11;}if(t)top=f;else if(b)bottom=f;else for(int j=0;j<6;j++){boolean s=true;for(int z:ids)s&=Math.abs(v[0][z-1]*Math.cos(j*Math.PI/3)+v[1][z-1]*Math.sin(j*Math.PI/3)-ap)<1e-11;if(s){sides[j]=f;break;}}}
                for(int j=0;j<v[0].length;j++)if(Math.abs(v[0][j]-ap)<1e-11&&Math.abs(v[1][j]-ap/Math.sqrt(3))<1e-11&&Math.abs(v[2][j]-H)<1e-11)ref=j+1;
                if(top<0||bottom<0||ref<0)throw new IllegalStateException("Flat hex face/reference invalid");for(int f:sides)if(f==0)throw new IllegalStateException("Missing periodic side");
                m.component("comp1").material().create("air","Common");m.component("comp1").material("air").propertyGroup("def").set("relpermittivity","1");m.component("comp1").material("air").propertyGroup("def").set("relpermeability","1");m.component("comp1").material("air").propertyGroup("def").set("electricconductivity","0[S/m]");
                m.component("comp1").physics().create("emw","ElectromagneticWaves","geom1");m.component("comp1").physics("emw").prop("ShapeProperty").set("order_electricfield","2");
                m.component("comp1").physics("emw").create("port1","Port",2);PhysicsFeature port=m.component("comp1").physics("emw").feature("port1");port.selection().set(new int[]{top});port.set("PortType","Periodic");port.set("PortName","1");port.set("PortExcitation","on");port.set("Pin","1[W]");
                double[]out=FlatPrecheckEntry.globalE(pol,theta,phi,true);String[]pout={Double.toString(out[0])+"[V/m]",Double.toString(out[1])+"[V/m]",Double.toString(out[2])+"[V/m]"};port.set("Polarization","UserDefined");port.set("InputType","E");port.set("E0",pout);port.set("Eampl",pout);port.set("alpha1_inc",theta+"[deg]");port.set("alpha2_inc",(phi+FlatPrecheckEntry.PORT_AZIMUTH_OFFSET_DEG)+"[deg]");
                port.create("ref1","PeriodicPortReferencePoint",0);port.feature("ref1").selection().set(new int[]{ref});
                for(int j=0;j<3;j++){String tag="pc"+j;m.component("comp1").physics("emw").create(tag,"PeriodicCondition",2);PhysicsFeature p=m.component("comp1").physics("emw").feature(tag);p.selection().set(new int[]{sides[j],sides[j+3]});p.set("PeriodicType","Floquet");p.set("Floquet_source","FromPeriodicPort");}
                m.component("comp1").physics("emw").create("imp","Impedance",2);m.component("comp1").physics("emw").feature("imp").selection().set(new int[]{bottom});m.component("comp1").physics("emw").feature("imp").set("surfaceImpedanceDefinition","userdefined");m.component("comp1").physics("emw").feature("imp").set("surfaceImpedance","Zs");
                m.study().create("std1");m.study("std1").create("freq","Frequency");m.study("std1").feature("freq").set("plist","f0");port.set("MaximumFrequency","UserDefined");port.set("fmax","f0");port.runCommand("addDiffractionOrders");
                List<String> others=new ArrayList<>();for(String tag:port.feature().tags()){PhysicsFeature f=port.feature(tag);if(f.getType().equals("OrthogonalPolarization"))others.add(f.getString("PortName"));else if(f.getType().equals("DiffractionOrder"))throw new IllegalStateException("Subwavelength control unexpectedly has nonzero open order");}
                if(others.size()!=1)throw new IllegalStateException("Flat control must include both zero-order polarizations");
                m.component("comp1").mesh().create("mesh1");MeshSequence ms=m.component("comp1").mesh("mesh1");ms.feature("size").set("custom","on");ms.feature("size").set("hmax",.35e-6);ms.feature("size").set("hmin",.1e-6);
                for(int j=0;j<3;j++){String tag="id"+j;ms.create(tag,"IdenticalMesh");ms.feature(tag).selection("group1").set(new int[]{sides[j]});ms.feature(tag).selection("group2").set(new int[]{sides[j+3]});}ms.create("ftet","FreeTet");ms.run();
                m.study("std1").createAutoSequences("all");String solver=null,variables=null;
                for(String s:m.sol().tags())for(String t:m.sol(s).feature().tags())if(m.sol(s).feature(t).getType().equals("Variables")){if(solver!=null)throw new IllegalStateException("Ambiguous flat variables solver");solver=s;variables=t;}
                if(solver==null)throw new IllegalStateException("Missing flat variables sequence");m.sol(solver).run(variables);
                double[] presolveK=new double[3],requested=FlatPrecheckEntry.requestedIncidentK(lambda,theta,phi);boolean preOK=true;
                String[] abc={"x","y","z"};for(int j=0;j<3;j++){presolveK[j]=glob(m,"real(emw.kInc"+abc[j]+"_1)");preOK&=Math.abs(presolveK[j]-requested[j])<.1;}
                System.out.println("CODEX_PRESOLVE_INCIDENT_K_JSON="+flatJson(flatObj("case_id",id,"actual_k",presolveK,"requested_k",requested,"passed",preOK)));System.out.flush();
                if(!preOK)throw new IllegalStateException("Flat compiled presolve k differs from requested global incidence");
                long start=System.nanoTime();m.sol(solver).runAll();double seconds=(System.nanoTime()-start)*1e-9;m.save(root+"\\"+id+"_solved_unevaluated.mph");
                double[] k=new double[3],want=FlatPrecheckEntry.requestedIncidentK(lambda,theta,phi);String[] axes={"x","y","z"};boolean wave=true;for(int j=0;j<3;j++){k[j]=integ(m,top,"real(emw.kInc"+axes[j]+"_1)")/area;wave&=Math.abs(k[j]-want[j])<.1;}
                Map<String,Object> fit=FlatPrecheckEntry.internalElectricPowerFit(m,area,0,H,lambda,k);
                double pin=((Number)fit.get("Pin_W")).doubleValue(),R=glob(m,"abs(emw.S11)^2");for(String pn:others){String compact="emw.S"+pn+"1",sep="emw.S"+pn+"_1";String expr=glob(m,"isdefined("+sep+")")>.5?sep:compact;R+=glob(m,"abs("+expr+")^2");}
                String h2="abs(emw.Hx)^2+abs(emw.Hy)^2+abs(emw.Hz)^2-abs(nx*emw.Hx+ny*emw.Hy+nz*emw.Hz)^2",e2="abs(emw.Ex)^2+abs(emw.Ey)^2+abs(emw.Ez)^2-abs(nx*emw.Ex+ny*emw.Ey+nz*emw.Ez)^2";
                double ah=integ(m,bottom,".5*real(Zs)*("+h2+")"),ae=integ(m,bottom,".5*real(1/Zs)*("+e2+")"),an=-integ(m,top,"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz"),ribc=glob(m,"R_ibc"),modeTotal=0,leak=0,side=0;
                for(int j=0;j<3;j++){double q=integ(m,top,"abs(emw.tEmode"+axes[j]+"_1)^2");modeTotal+=q;}
                double cp=Math.cos(Math.toRadians(phi)),sp=Math.sin(Math.toRadians(phi)),modeAxisX=pol.equals("s")?-sp:cp,modeAxisY=pol.equals("s")?cp:sp;leak=1-integ(m,top,"abs(("+modeAxisX+")*emw.tEmodex_1+("+modeAxisY+")*emw.tEmodey_1)^2")/modeTotal;
                for(int j=0;j<3;j++)side=Math.max(side,Math.abs(integ(m,sides[j],"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz")+integ(m,sides[j+3],"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz")));
                boolean fitOK=Boolean.TRUE.equals(fit.get("condition_passed"))&&Boolean.TRUE.equals(fit.get("residual_passed"))&&Boolean.TRUE.equals(fit.get("LOO_Pin_passed"))&&Boolean.TRUE.equals(fit.get("all_cut_areas_match_geometry"));
                Map<String,Object> checks=flatObj("requested_incident_wavevector",wave,"source_mode_purity",Double.isFinite(modeTotal)&&modeTotal>0&&Double.isFinite(leak)&&leak>=-1e-12&&leak<=1e-6,"independent_incident_power",Math.abs(pin-1)<=.005,"shared_source_fit",fitOK,"analytic_IBC_R",Math.abs(R-ribc)<=.002,"port_vs_Ti_loss_absolute",Math.abs(1-R-ah)<=.005,"port_vs_Ti_loss_relative",Math.abs(1-R-ah)/Math.max(1-R,1e-12)<=.02,"Ti_H_vs_E_relative",Math.abs(ah-ae)/Math.max(ah,1e-12)<=.02,"net_vs_Ti_loss_absolute",Math.abs(an-ah)<=.005,"periodic_sides",side<=.002);
                checks.put("Ti_H_vs_E_absolute",Math.abs(ah-ae)<=.005);checks.put("finite",Double.isFinite(pin)&&Double.isFinite(R)&&Double.isFinite(ah)&&Double.isFinite(ae)&&Double.isFinite(an)&&Double.isFinite(ribc)&&Double.isFinite(leak)&&Double.isFinite(side));
                boolean pass=true;for(Object b:checks.values())pass&=Boolean.TRUE.equals(b);m.save(root+"\\"+id+"_model.mph");
                Map<String,Object> row=flatObj("case_id",id,"state",pass?"PASS":"FAILED_FLAT_ADMISSION","wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"polarization",pol,"solver_executed",true,"source_fit_passed",fitOK,"incident_power_passed",Math.abs(pin-1)<=.005,"source_wavevector_and_polarization_passed",wave&&leak<=1e-6,"analytic_IBC_R_passed",Math.abs(R-ribc)<=.002,"actual_R",R,"analytic_IBC_R",ribc,"actual_Pin_W",pin,"A_Ti_H",ah,"A_Ti_E",ae,"A_net",an,"actual_incident_k",k,"requested_incident_k",want,"mode_leakage",leak,"side_pair_mismatch",side,"checks",checks,"shared_production_power_fit",fit,"model_filename",id+"_model.mph","solve_seconds",seconds,"purpose","SOURCE_NORMALIZATION_IMPLEMENTATION_ONLY_NOT_CAVITY_PERFORMANCE","finished_utc",Instant.now().toString());
                row.put("solver_completed",true);row.put("requested_global_source_E_out",Arrays.asList(pout));row.put("source_mode_target_tangential_component",pol.equals("s")?"globaltangentials=(-sinphi,cosphi,0)":"globaltangentialp=(cosphi,sinphi,0)");row.put("normalized_source_mode_target_tangential_unit_axis",new double[]{modeAxisX,modeAxisY,0});row.put("mode_total",modeTotal);row.put("R_total",R);row.put("internal_E_power_fit",fit);row.put("presolve_k_verified",preOK);
                System.out.println("CODEX_FLAT_ADMISSION_JSON="+flatJson(row));System.out.flush();if(!pass)throw new IllegalStateException("Native flat source/power/IBC control failed: "+id);
            }finally{if(m!=null)ModelUtil.remove("FlatAdmission");}
        }
        System.out.println("CODEX_ALL_FLAT_CONTROLS_PASSED=true");
    }

public static void main(String[]args)throws Exception{throw new IllegalStateException("Compile-only;neverexecute");}
}
