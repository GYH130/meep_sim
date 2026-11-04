import com.comsol.model.*;
import com.comsol.model.physics.PhysicsFeature;
import com.comsol.model.util.ModelUtil;
import java.time.Instant;
import java.util.*;

/** Tiny planar Ti-IBC controls for the shared production source/power fit.
 * Subwavelength hexagon has the same reference-vertex convention, not the
 * research cavity geometry. These results admit implementation only.
 */
public class FlatObliqueAdmission {
    private static int ev=0;
    private static Map<String,Object> obj(Object... kv){Map<String,Object> r=new LinkedHashMap<>();for(int i=0;i<kv.length;i+=2)r.put((String)kv[i],kv[i+1]);return r;}
    private static String json(Object x){
        if(x==null)return "null";if(x instanceof String)return "\""+((String)x).replace("\\","\\\\").replace("\"","\\\"").replace("\n","\\n").replace("\r","\\r")+"\"";
        if(x instanceof Number){if(!Double.isFinite(((Number)x).doubleValue()))throw new IllegalStateException("Nonfinite calibration metric");return x.toString();}
        if(x instanceof Boolean)return x.toString();StringJoiner j;
        if(x instanceof Map){j=new StringJoiner(",","{","}");for(Object k:((Map<?,?>)x).keySet())j.add(json(k.toString())+":"+json(((Map<?,?>)x).get(k)));return j.toString();}
        if(x instanceof Collection){j=new StringJoiner(",","[","]");for(Object q:(Collection<?>)x)j.add(json(q));return j.toString();}
        if(x instanceof double[]){List<Double> r=new ArrayList<>();for(double q:(double[])x)r.add(q);return json(r);}throw new IllegalArgumentException("Unsupported calibration JSON");
    }
    private static double integ(Model m,int face,String e){String tag="flatEv"+(++ev);m.result().numerical().create(tag,"IntSurface");try{NumericalFeature n=m.result().numerical(tag);n.selection().set(new int[]{face});n.set("data","dset1");n.set("expr",e);n.set("intorderactive","on");n.set("intorder",8);return n.getReal()[0][0];}finally{m.result().numerical().remove(tag);}}
    private static double glob(Model m,String e){String tag="flatEv"+(++ev);m.result().numerical().create(tag,"EvalGlobal");try{NumericalFeature n=m.result().numerical(tag);n.set("data","dset1");n.set("expr",e);return n.getReal()[0][0];}finally{m.result().numerical().remove(tag);}}
    public static void runPoint(String root,String fixedId,double lambdaIn,double thetaIn,double ntiIn,double ktiIn,String pol,double phi) throws Exception {
        System.out.println("CODEX_ARITHMETIC_SELF_TEST_JSON="+json(CavityFloorFilm.pureArithmeticSelfTest()));System.out.flush();
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
            String id=fixedId;Model m=null;ev=0;CavityFloorFilm.setProbeContext(id);
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
                double[]out=CavityFloorFilm.globalE(pol,theta,phi,true);String[]pout={Double.toString(out[0])+"[V/m]",Double.toString(out[1])+"[V/m]",Double.toString(out[2])+"[V/m]"};port.set("Polarization","UserDefined");port.set("InputType","E");port.set("E0",pout);port.set("Eampl",pout);port.set("alpha1_inc",theta+"[deg]");port.set("alpha2_inc",(phi+CavityFloorFilm.PORT_AZIMUTH_OFFSET_DEG)+"[deg]");
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
                double[] presolveK=new double[3],requested=CavityFloorFilm.requestedIncidentK(lambda,theta,phi);boolean preOK=true;
                String[] abc={"x","y","z"};for(int j=0;j<3;j++){presolveK[j]=glob(m,"real(emw.kInc"+abc[j]+"_1)");preOK&=Math.abs(presolveK[j]-requested[j])<.1;}
                System.out.println("CODEX_PRESOLVE_INCIDENT_K_JSON="+json(obj("case_id",id,"actual_k",presolveK,"requested_k",requested,"passed",preOK)));System.out.flush();
                if(!preOK)throw new IllegalStateException("Flat compiled presolve k differs from requested global incidence");
                long start=System.nanoTime();m.sol(solver).runAll();double seconds=(System.nanoTime()-start)*1e-9;m.save(root+"\\"+id+"_solved_unevaluated.mph");
                double[] k=new double[3],want=CavityFloorFilm.requestedIncidentK(lambda,theta,phi);String[] axes={"x","y","z"};boolean wave=true;for(int j=0;j<3;j++){k[j]=integ(m,top,"real(emw.kInc"+axes[j]+"_1)")/area;wave&=Math.abs(k[j]-want[j])<.1;}
                Map<String,Object> fit=CavityFloorFilm.internalElectricPowerFit(m,area,0,H,lambda,k);
                double pin=((Number)fit.get("Pin_W")).doubleValue(),R=glob(m,"abs(emw.S11)^2");for(String pn:others){String compact="emw.S"+pn+"1",sep="emw.S"+pn+"_1";String expr=glob(m,"isdefined("+sep+")")>.5?sep:compact;R+=glob(m,"abs("+expr+")^2");}
                String h2="abs(emw.Hx)^2+abs(emw.Hy)^2+abs(emw.Hz)^2-abs(nx*emw.Hx+ny*emw.Hy+nz*emw.Hz)^2",e2="abs(emw.Ex)^2+abs(emw.Ey)^2+abs(emw.Ez)^2-abs(nx*emw.Ex+ny*emw.Ey+nz*emw.Ez)^2";
                double ah=integ(m,bottom,".5*real(Zs)*("+h2+")"),ae=integ(m,bottom,".5*real(1/Zs)*("+e2+")"),an=-integ(m,top,"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz"),ribc=glob(m,"R_ibc"),modeTotal=0,leak=0,side=0;
                for(int j=0;j<3;j++){double q=integ(m,top,"abs(emw.tEmode"+axes[j]+"_1)^2");modeTotal+=q;}
                double cp=Math.cos(Math.toRadians(phi)),sp=Math.sin(Math.toRadians(phi)),modeAxisX=pol.equals("s")?-sp:cp,modeAxisY=pol.equals("s")?cp:sp;leak=1-integ(m,top,"abs(("+modeAxisX+")*emw.tEmodex_1+("+modeAxisY+")*emw.tEmodey_1)^2")/modeTotal;
                for(int j=0;j<3;j++)side=Math.max(side,Math.abs(integ(m,sides[j],"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz")+integ(m,sides[j+3],"emw.Poavx*nx+emw.Poavy*ny+emw.Poavz*nz")));
                boolean fitOK=Boolean.TRUE.equals(fit.get("condition_passed"))&&Boolean.TRUE.equals(fit.get("residual_passed"))&&Boolean.TRUE.equals(fit.get("LOO_Pin_passed"))&&Boolean.TRUE.equals(fit.get("all_cut_areas_match_geometry"));
                Map<String,Object> checks=obj("requested_incident_wavevector",wave,"source_mode_purity",Double.isFinite(modeTotal)&&modeTotal>0&&Double.isFinite(leak)&&leak>=-1e-12&&leak<=1e-6,"independent_incident_power",Math.abs(pin-1)<=.005,"shared_source_fit",fitOK,"analytic_IBC_R",Math.abs(R-ribc)<=.002,"port_vs_Ti_loss_absolute",Math.abs(1-R-ah)<=.005,"port_vs_Ti_loss_relative",Math.abs(1-R-ah)/Math.max(1-R,1e-12)<=.02,"Ti_H_vs_E_relative",Math.abs(ah-ae)/Math.max(ah,1e-12)<=.02,"net_vs_Ti_loss_absolute",Math.abs(an-ah)<=.005,"periodic_sides",side<=.002);
                checks.put("Ti_H_vs_E_absolute",Math.abs(ah-ae)<=.005);checks.put("finite",Double.isFinite(pin)&&Double.isFinite(R)&&Double.isFinite(ah)&&Double.isFinite(ae)&&Double.isFinite(an)&&Double.isFinite(ribc)&&Double.isFinite(leak)&&Double.isFinite(side));
                boolean pass=true;for(Object b:checks.values())pass&=Boolean.TRUE.equals(b);m.save(root+"\\"+id+"_model.mph");
                Map<String,Object> row=obj("case_id",id,"state",pass?"PASS":"FAILED_FLAT_ADMISSION","wavelength_um",lambda,"theta_deg",theta,"phi_global_deg",phi,"polarization",pol,"solver_executed",true,"source_fit_passed",fitOK,"incident_power_passed",Math.abs(pin-1)<=.005,"source_wavevector_and_polarization_passed",wave&&leak<=1e-6,"analytic_IBC_R_passed",Math.abs(R-ribc)<=.002,"actual_R",R,"analytic_IBC_R",ribc,"actual_Pin_W",pin,"A_Ti_H",ah,"A_Ti_E",ae,"A_net",an,"actual_incident_k",k,"requested_incident_k",want,"mode_leakage",leak,"side_pair_mismatch",side,"checks",checks,"shared_production_power_fit",fit,"model_filename",id+"_model.mph","solve_seconds",seconds,"purpose","SOURCE_NORMALIZATION_IMPLEMENTATION_ONLY_NOT_CAVITY_PERFORMANCE","finished_utc",Instant.now().toString());
                row.put("solver_completed",true);row.put("requested_global_source_E_out",Arrays.asList(pout));row.put("source_mode_target_tangential_component",pol.equals("s")?"globaltangentials=(-sinphi,cosphi,0)":"globaltangentialp=(cosphi,sinphi,0)");row.put("normalized_source_mode_target_tangential_unit_axis",new double[]{modeAxisX,modeAxisY,0});row.put("mode_total",modeTotal);row.put("R_total",R);row.put("internal_E_power_fit",fit);row.put("presolve_k_verified",preOK);
                System.out.println("CODEX_FLAT_ADMISSION_JSON="+json(row));System.out.flush();if(!pass)throw new IllegalStateException("Native flat source/power/IBC control failed: "+id);
            }finally{if(m!=null)ModelUtil.remove("FlatAdmission");}
        }
        System.out.println("CODEX_ALL_FLAT_CONTROLS_PASSED=true");
    }
}
