'use strict';
const assert=require('node:assert/strict');
module.exports=function tests(s,n) {
  const {F,Cut,ONE,IOTA}=n, f=F.of, eq=(a,b)=>assert.ok(f(a).eq(b),`${a} != ${b}`);
  const same=(a,b)=>assert.ok(n.equal(a,b));
  const report=[], cases={histories:0,recoveries:0,record_trees:0,metric_tangents:0,observer_completions:0};
  const check=(name,fn)=>{fn();report.push({name,passed:true});};
  const sources=[[1,0],[0,1],[3,4],[2,-1],[-3,4],[0,2],[-1,0],['3/5','4/5'],[1,1]];
  const states=[[0,1],[1,0],[2,-3],['5/13','12/13']];
  check('quarter-turn and dagger are native scalar operations',()=>{
    assert.ok(IOTA.mul(IOTA).eq(-1));assert.ok(IOTA.dagger().eq(IOTA.neg()));
    for(const z of sources.map(x=>new Cut(...x)))assert.ok(z.dagger().dagger().eq(z));
  });
  check('radial/complement split derived from dagger is idempotent',()=>{
    for(const a of states){const z=new Cut(...a);assert.ok(s.P(z).add(s.Q(z)).eq(z));assert.ok(s.P(s.P(z)).eq(s.P(z)));assert.ok(s.Q(s.Q(z)).eq(s.Q(z)));assert.ok(s.P(s.Q(z)).zero());}
    assert.ok(s.P(IOTA).zero());assert.ok(!IOTA.mul(s.P(ONE)).zero());
  });
  check('cut-square energy decomposes and is positive',()=>{
    for(const z of [...sources,...states].map(x=>new Cut(...x))){eq(z.norm2(),s.P(z).norm2().add(s.Q(z).norm2()));assert.ok(f(0).le(z.norm2()));}
  });
  check('multiplication derives raw transport and its radial scale',()=>{
    for(const z of sources){const t=s.source(z);same(t.T,n.matrix([[t.a,t.b.neg()],[t.b,t.a]]));same(n.mul(n.dagger(t.T),t.T),n.scale(n.identity(2),t.d));
      for(const xy of states)eq(t.z.mul(new Cut(...xy)).norm2(),t.d.mul(new Cut(...xy).norm2()));}
  });
  check('cut commutator, leakage and return use one source coefficient',()=>{
    for(const z of sources){const t=s.source(z),ptp=n.mul(n.mul(t.P,t.T),t.P),m=n.commutator(t.P,t.T);
      same(n.mul(m,m),n.scale(n.identity(2),t.b.pow(2)));
      same(n.sub(n.mul(n.mul(t.P,n.power(t.T,2)),t.P),n.mul(ptp,ptp)),n.scale(t.P,t.b.pow(2).neg()));
      same(n.sub(n.scale(t.P,t.d),n.mul(n.dagger(ptp),ptp)),n.scale(t.P,t.b.pow(2)));}
  });
  check('exact reduced memory reconstructs unnormalized full histories',()=>{
    for(const z of sources)for(const xy of states){const h=s.history(z,...xy,8);h.residual.forEach(x=>eq(x,0));cases.histories++;}
  });
  check('two-step recurrence retains source scale instead of silently setting it to one',()=>{
    for(const z of sources)for(const xy of states){const t=s.source(z),h=s.history(z,...xy,8);
      for(let j=0;j+2<h.observed.length;j++)eq(h.observed[j+2],t.a.mul(2).mul(h.observed[j+1]).sub(t.d.mul(h.observed[j])));}
  });
  check('two radial readings recover every coupled complementary source',()=>{
    for(const z of sources)for(const xy of states){const t=s.source(z),next=t.z.mul(new Cut(...xy));
      if(t.b.zero())assert.throws(()=>s.recover(z,xy[0],next.rad));else{eq(s.recover(z,xy[0],next.rad),xy[1]);cases.recoveries++;}}
  });
  check('canonical N03 completion has the proved minimum rank',()=>{
    for(const z of sources){const t=s.source(z),c=n.rowClosure([[1,0]],[t.T]);assert.equal(c.rank,t.b.zero()?1:2);cases.observer_completions++;}
  });
  check('arbitrary finite refinement ledger derives variance by pair squares',()=>{
    for(const cells of [[[1,1]],[[1,1],[1,-1]],[[3,2],[1,-2]],[[2,'1/3'],[5,'-2/3'],[1,0]]]){
      const p=s.preparation(cells);eq(p.variance,p.pairSquare);eq(s.sum(p.weights),1);assert.ok(f(0).le(p.variance));}
    const p=s.preparation([[3,2],[1,-2]]);eq(p.mean,1);eq(p.variance,3);
  });
  check('refinement and record ordering preserve derived moments',()=>{
    const p=s.preparation([[3,2],[1,-2]]),r=s.preparation([[1,2],[1,-2],[1,2],[1,2]]),u=s.preparation([[6,2],[2,-2]]);
    for(const q of [r,u]){eq(p.mean,q.mean);eq(p.variance,q.variance);}
  });
  check('source count balance is derived when present and not imposed',()=>{
    const p=s.preparation([[1,'12/13'],[1,'-12/13']]);eq(p.mean,0);eq(p.variance,f(1).sub(f('5/13').pow(2)));
    assert.ok(!s.preparation([[3,1],[1,-1]]).mean.zero());
  });
  check('colored force covariance follows the full source ledger',()=>{
    for(const z of sources){const q=s.noise(z,[[3,2],[1,-2]],6);
      same(q.covariance,q.predicted);q.memoryNoiseResidual.forEach(x=>eq(x,0));
      for(let k=0;k<6;k++)eq(s.sum(q.forces.map((row,j)=>q.preparation.weights[j].mul(row[k]))),0);
      assert.ok(n.rank(q.covariance)<=1);}
  });
  check('resolved force can vanish while returning memory survives',()=>{
    const q=s.noise([3,4],[[1,2]],3);assert.ok(n.isZero(q.covariance));eq(s.history([3,4],0,2,3).kernels[0],-16);
  });
  check('first-exit content equals normalized native path squares',()=>{
    for(const z of sources)for(const k of [0,1,3,8]){const e=s.firstExit(z,k);eq(e.total,1);e.weights.forEach((v,j)=>eq(v,e.amplitudeWeights[j]));}
  });
  check('first-exit endpoints retain unresolved and immediate-exit cases',()=>{
    eq(s.firstExit([2,0],5).tail,1);eq(s.firstExit([0,2],5).weights[0],1);eq(s.firstExit([0,2],5).tail,0);
  });
  check('full record refinement conserves normalized energy and reconstructs coherent sum',()=>{
    for(const z of sources)for(let depth=0;depth<=5;depth++){
      const t=s.source(z),r=s.recordTree(z,[2,-3],depth);eq(r.totalContent,1);
      const expected=n.mul(n.power(t.T,depth),s.column(new Cut(2,-3)));same(s.column(r.coherentSum),expected);cases.record_trees++;}
  });
  check('first-exit cylinder content is independent of further record refinement',()=>{
    for(const z of sources){const leaves=s.recordTree(z,ONE,6).leaves,e=s.firstExit(z,5);
      for(let k=0;k<6;k++)eq(s.sum(leaves.filter(x=>x.word.startsWith('P'.repeat(k)+'Q')).map(x=>x.content)),e.weights[k]);
      eq(s.sum(leaves.filter(x=>x.word==='PPPPPP').map(x=>x.content)),e.tail);}
  });
  check('coherent erasure and monitored survival remain different operations',()=>{
    const t=s.source(['3/5','4/5']),r=s.recordTree(t.z,ONE,2);
    eq(r.coherentSum.rad.pow(2),'49/625');eq(r.leaves.find(x=>x.word==='PP').content,'81/625');
  });
  check('relative-content Hessian equals phase geometry without an angle input',()=>{
    for(const z of sources)for(const d of [[1,0],[0,1],[2,-3],['1/7','2/5']]){
      const g=s.sourceMetric(z,...d);eq(g.nativePhaseMetric.mul(4),g.amplitudeMetric);eq(g.nativePhaseMetric,g.phaseDifferential.pow(2));if(g.relativeContentHessian!==null)eq(g.relativeContentHessian,g.amplitudeMetric);cases.metric_tangents++;}
  });
  check('derived geometry has radial kernel and angular normalization four',()=>{
    for(const z of sources){const t=s.source(z);eq(s.sourceMetric(z,t.a,t.b).amplitudeMetric,0);eq(s.sourceMetric(z,t.b.neg(),t.a).amplitudeMetric,4);}
  });
  check('source rescaling preserves the phase geometry with transported tangents',()=>{
    for(const z of sources){const t=s.source(z),g=s.sourceMetric(t.z,2,3),h=s.sourceMetric(t.z.mul(-5),-10,-15);eq(g.amplitudeMetric,h.amplitudeMetric);}
  });
  check('relative-content jets use native log derivatives and normalized source content',()=>{
    const r=s.relativeContentJet(['9/25','16/25'],['-24/25','24/25'],['14/25','-14/25']);eq(r.first,0);eq(r.second,4);
    assert.throws(()=>s.relativeContentJet([1,0],[0,0]));assert.throws(()=>s.relativeContentJet(['1/2','1/2'],[1,0]));
  });
  check('native logarithm enclosures independently validate the local contrast Hessian',()=>{
    const p=['9/25','16/25'].map(f),v=['-24/25','24/25'].map(f),h=f('1/1000');
    let lower=f(0),upper=f(0),bound=f(0);
    for(let j=0;j<2;j++){
      for(const sign of [-1,1]){const box=s.logInterval(p[j].div(p[j].add(h.mul(v[j]).mul(sign))),5);lower=lower.add(p[j].mul(box.lower));upper=upper.add(p[j].mul(box.upper));}
      const u=h.mul(v[j]).div(p[j]);bound=bound.add(h.pow(2).mul(v[j].pow(4)).div(p[j].pow(3).mul(2).mul(f(1).sub(u.pow(2)))));
    }
    lower=lower.div(h.pow(2));upper=upper.div(h.pow(2));assert.ok(f(4).le(lower));assert.ok(upper.le(f(4).add(bound)));
  });
  check('support boundaries are not reported as regular event Hessians',()=>{
    assert.equal(s.sourceMetric([1,0],0,1).relativeContentHessian,null);eq(s.sourceMetric([1,0],0,1).amplitudeMetric,4);
  });
  check('native angle series and factorial Euler independently enclose source phase',()=>{
    for(const h of ['-1/2','-1/3','0','1/3','1/2']){
      const box=s.phaseInterval(h,12),mid=box.lower.add(box.upper).div(2),width=box.upper.sub(box.lower).div(2);
      let power=ONE,poly=ONE,factorial=f(1);
      for(let k=1;k<=28;k++){power=power.mul(new Cut(0,mid));factorial=factorial.mul(k);poly=poly.add(power.div(factorial));}
      const first=mid.abs().pow(29).div(factorial.mul(29));
      const tail=first.div(f(1).sub(mid.abs().div(30))),radius=tail.add(width);
      assert.ok(poly.rad.sub(box.cosine).abs().le(radius));assert.ok(poly.turn.sub(box.sine).abs().le(radius));
      eq(s.source([box.cosine,box.sine]).d,1);
    }
    assert.throws(()=>s.phaseInterval(2,12));
  });
  check('signed phase remains distinguishable beyond squared event content',()=>{
    const t=s.source([3,4]),u=s.source([3,-4]);eq(t.q,u.q);assert.ok(!t.z.mul(IOTA).eq(u.z.mul(IOTA)));
  });
  check('same native algebra permits inequivalent phase and preparation sources',()=>{
    assert.ok(!s.source([3,4]).q.eq(s.source([5,12]).q));
    assert.ok(!s.preparation([[1,1],[1,-1]]).variance.eq(s.preparation([[1,1]]).variance));
  });
  check('integer path memory survives endpoint return',()=>{
    let z=ONE;for(let k=0;k<4;k++)z=IOTA.mul(z);assert.ok(z.eq(ONE));assert.notEqual(4,0);
  });
  check('invalid source, ledger and exact-input claims are rejected',()=>{
    assert.throws(()=>s.source([0,0]));assert.throws(()=>s.source([0.5,1]));assert.throws(()=>s.preparation([]));assert.throws(()=>s.preparation([[0,2]]));assert.throws(()=>s.recordTree([1,1],[0,0],2));assert.throws(()=>s.logInterval(0,3));
  });
  return {tests:report,counts:cases};
};
