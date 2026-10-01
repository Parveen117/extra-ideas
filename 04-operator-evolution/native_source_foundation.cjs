'use strict';
/** R14 adapter: all arithmetic is supplied by the pinned canonical RKF engine.
 * The input is a nonzero native source scalar plus a finite preparation ledger.
 * No transcendental library, random generator or physical measurement axiom.
 * Finite exact witnesses do not constitute universal proof or event actualization.
 */
module.exports = function foundation(n) {
  const {F, Cut, ONE, IOTA} = n;
  const f = x => F.of(x);
  const sum = xs => xs.reduce((a,b) => a.add(b), f(0));
  const cut = x => Cut.of(x);
  function count(k) {
    if (!Number.isSafeInteger(k) || k < 0) throw new RangeError('Nonnegative finite count required');
    return k;
  }
  const P = z => cut(z).add(cut(z).dagger()).div(2);
  const Q = z => cut(z).sub(cut(z).dagger()).div(2);
  const column = z => [[cut(z).rad],[cut(z).turn]];
  const coords = op => {
    const a = op(ONE), b = op(IOTA);
    return n.matrix([[a.rad,b.rad],[a.turn,b.turn]]);
  };
  function source(z) {
    z = cut(z);
    if (z.zero()) throw new RangeError('Zero source is a separate stratum; phase and normalized content undefined');
    const a=z.rad, b=z.turn, d=z.norm2();
    return {z,a,b,d, r:a.pow(2).div(d), q:b.pow(2).div(d),
      J:coords(v=>v.dagger()), P:coords(P), Q:coords(Q),
      R:coords(v=>IOTA.mul(v)), T:coords(v=>z.mul(v))};
  }
  function history(z,x,y,length) {
    const t=source(z); x=f(x);y=f(y);count(length);
    let v=new Cut(x,y); const full=[v], observed=[x], predicted=[], residual=[];
    for(let k=0;k<length;k++) {
      const memory=sum(observed.slice(0,k).map((u,j)=>t.b.pow(2).neg().mul(t.a.pow(k-1-j)).mul(u)));
      const force=t.b.neg().mul(t.a.pow(k)).mul(y);
      const estimate=t.a.mul(observed[k]).add(memory).add(force);
      v=t.z.mul(v);full.push(v);observed.push(v.rad);predicted.push(estimate);residual.push(v.rad.sub(estimate));
    }
    return {full,observed,predicted,residual,
      kernels:Array.from({length},(_,k)=>t.b.pow(2).neg().mul(t.a.pow(k)))};
  }
  function preparation(cells) {
    if(!Array.isArray(cells)||cells.length===0)throw new RangeError('Nonempty source refinement ledger required');
    const rows=cells.map(([multiplicity,y])=>{
      if(!Number.isSafeInteger(multiplicity)||multiplicity<1)throw new RangeError('Positive integer cell multiplicity required');
      return {count:f(multiplicity),y:f(y)};
    });
    const total=sum(rows.map(x=>x.count));
    const weights=rows.map(x=>x.count.div(total));
    const mean=sum(rows.map((x,j)=>weights[j].mul(x.y)));
    const variance=sum(rows.map((x,j)=>weights[j].mul(x.y.sub(mean).pow(2))));
    const pairSquare=sum(rows.flatMap((x,j)=>rows.map((y,k)=>weights[j].mul(weights[k]).mul(x.y.sub(y.y).pow(2))))).div(2);
    return {rows,total,weights,mean,variance,pairSquare};
  }
  function noise(z,cells,length) {
    const t=source(z), s=preparation(cells);count(length);
    const forces=s.rows.map(x=>Array.from({length},(_,k)=>t.b.neg().mul(t.a.pow(k)).mul(x.y.sub(s.mean))));
    const covariance=Array.from({length},(_,k)=>Array.from({length},(_,l)=>sum(forces.map((row,j)=>s.weights[j].mul(row[k]).mul(row[l])))));
    const predicted=Array.from({length},(_,k)=>Array.from({length},(_,l)=>t.b.pow(2).mul(t.a.pow(k+l)).mul(s.variance)));
    return {preparation:s,forces,covariance,predicted,
      memoryNoiseResidual:Array.from({length},(_,k)=>covariance[k][0].sub(t.b.pow(2).mul(t.a.pow(k)).mul(s.variance)))};
  }
  function recover(z,x0,x1) {
    const t=source(z);if(t.b.zero())throw new RangeError('Complement is invisible to this forward observer');
    return t.a.mul(x0).sub(x1).div(t.b);
  }
  function firstExit(z,lastIndex) {
    const t=source(z);count(lastIndex);
    const weights=Array.from({length:lastIndex+1},(_,k)=>t.q.mul(t.r.pow(k)));
    const amplitudeWeights=Array.from({length:lastIndex+1},(_,k)=>t.b.pow(2).mul(t.a.pow(2*k)).div(t.d.pow(k+1)));
    const tail=t.r.pow(lastIndex+1);
    return {weights,amplitudeWeights,tail,total:sum(weights).add(tail),r:t.r,q:t.q};
  }
  function recordTree(z,initial,depth) {
    const t=source(z), v=cut(initial);count(depth);
    if(v.zero())throw new RangeError('Nonzero preparation required for normalized record content');
    if(depth>12)throw new RangeError('Finite audit tree budget exceeded');
    let leaves=[{word:'',amplitude:v}];
    for(let k=0;k<depth;k++)leaves=leaves.flatMap(row=>{
      const next=t.z.mul(row.amplitude);
      return [{word:row.word+'P',amplitude:P(next)},{word:row.word+'Q',amplitude:Q(next)}];
    });
    const totalEnergy=t.d.pow(depth).mul(v.norm2());
    leaves=leaves.map(x=>({...x,content:x.amplitude.norm2().div(totalEnergy)}));
    return {leaves,totalEnergy,totalContent:sum(leaves.map(x=>x.content)),
      coherentSum:leaves.reduce((a,b)=>a.add(b.amplitude),new Cut())};
  }
  function sourceMetric(z,da,db) {
    const t=source(z);da=f(da);db=f(db);
    const dn=t.a.mul(da).add(t.b.mul(db)).mul(2);
    const dr=t.a.mul(da).mul(2).mul(t.d).sub(t.a.pow(2).mul(dn)).div(t.d.pow(2));
    const dq=dr.neg();
    const phaseNumerator=t.a.mul(db).sub(t.b.mul(da));
    const amplitudeMetric=phaseNumerator.pow(2).mul(4).div(t.d.pow(2));
    const nativePhaseMetric=da.pow(2).add(db.pow(2)).div(t.d).sub(t.a.mul(da).add(t.b.mul(db)).pow(2).div(t.d.pow(2)));
    const regular=!t.r.zero()&&!t.q.zero();
    const relativeContentHessian=regular?dr.pow(2).div(t.r).add(dq.pow(2).div(t.q)):null;
    return {dr,dq,phaseDifferential:phaseNumerator.div(t.d),nativePhaseMetric,amplitudeMetric,relativeContentHessian,
      support:regular?'regular':'boundary: amplitude extension only'};
  }
  function relativeContentJet(probabilities,first,second=first.map(()=>0)) {
    const p=probabilities.map(f),v=first.map(f),w=second.map(f);
    if(p.length!==v.length||p.length!==w.length||!sum(p).eq(1)||!sum(v).zero()||!sum(w).zero()||p.some(x=>x.le(0)))
      throw new RangeError('Positive normalized content and normalized first/second jets required');
    return {first:sum(v).neg(),second:sum(p.map((x,j)=>v[j].pow(2).div(x).sub(w[j])))};
  }
  function logInterval(x,order) {
    x=f(x);count(order);if(x.le(0))throw new RangeError('Positive radial log input required');
    const u=x.sub(1).div(x.add(1));
    const center=sum(Array.from({length:order+1},(_,k)=>u.pow(2*k+1).mul(2).div(2*k+1)));
    const radius=u.abs().pow(2*order+3).mul(2).div(f(2*order+3).mul(f(1).sub(u.pow(2))));
    return {lower:center.sub(radius),upper:center.add(radius),center,radius};
  }
  function phaseInterval(h,order) {
    h=f(h);count(order);if(!h.abs().le(1))throw new RangeError('Use a native Cayley subchart with |h| <= 1');
    const partial=sum(Array.from({length:order+1},(_,k)=>h.pow(2*k+1).mul(k%2?-2:2).div(2*k+1)));
    const next=partial.add(h.pow(2*order+3).mul((order+1)%2?-2:2).div(2*order+3));
    return {lower:partial.le(next)?partial:next,upper:partial.le(next)?next:partial,
      cosine:f(1).sub(h.pow(2)).div(f(1).add(h.pow(2))),sine:h.mul(2).div(f(1).add(h.pow(2)))};
  }
  return {f,sum,P,Q,column,coords,source,history,preparation,noise,recover,firstExit,recordTree,sourceMetric,relativeContentJet,logInterval,phaseInterval};
};
