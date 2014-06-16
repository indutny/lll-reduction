var bn = require('bn.js');

function LLL() {
}
module.exports = LLL;

LLL.reduce = function reduce(basis) {
  return new LLL().reduce(basis);
};

LLL.prototype.vec = function vec(list) {
  if (Array.isArray(list))
    return new Vector(list);
  else
    return list;
};

LLL.prototype.basis = function basis(b) {
  var r = new Array(b.length);
  for (var i = 0; i < b.length; i++)
    r[i] = this.vec(b[i]);
  return r;
};

LLL.prototype.gramSchmidt = function gramSchmidt(basis) {
  var b = this.basis(basis);
  var n = basis[0].length;
  var res = new Array(b.length);
  var coeff = new Array(b.length);
  var squares = new Array(b.length);

  // First vector is kept unchanged
  res[0] = b[0];
  squares[0] = res[0].dot(res[0]);
  coeff[0] = [];
  for (var i = 1; i < b.length; i++) {
    var v = b[i];
    var co = new Array(i);

    var r = v;

    // Subtract projections from the vector
    for (var j = 0; j < i; j++) {
      var u = res[j];
      co[j] = u.dot(v).div(squares[j]);
      r = r.sub(u.mul(co[j]));
    }
    res[i] = r;
    coeff[i] = co;
    squares[i] = r.dot(r);
  }

  return {
    basis: res,
    coeff: coeff,
    squares: squares
  };
};

LLL.prototype.sizeReduce = function sizeReduce(i, j, b, coeff) {
  var q = new Fraction(coeff[i][j].toBN(true));
  if (q.cmpn(0) === 0)
    return;

  b[i] = b[i].sub(b[j].mul(q));
  for (var k = 0; k < j; k++)
    coeff[i][k] = coeff[i][k].sub(q.mul(coeff[j][k]));
  coeff[i][k] = coeff[i][k].sub(q);
};

LLL.prototype.reduce = function reduce(input, dn, dd) {
  var n = input.length;
  var d = input[0].length;
  var b = this.basis(input);

  var gram = this.gramSchmidt(b);
  var coeff = gram.coeff;
  var squares = gram.squares;

  var k = 1;
  while (k < n) {
    this.sizeReduce(k, k - 1, b, coeff);

    // Lovasz condition
    var lhs = squares[k];
    var rhs = new Fraction(dn, dd).sub(coeff[k][k - 1].sqr())
                                  .mul(squares[k - 1]);
    if (lhs.sub(rhs).cmpn(0) < 0) {
      // Swap (k - 1) and k
      var t = b[k];
      b[k] = b[k - 1];
      b[k - 1] = t;

      // Update Gram coefficients
      // TODO(indutny): use optimal update method
      var gram = this.gramSchmidt(b);
      var coeff = gram.coeff;
      var squares = gram.squares;

      k = Math.max(k - 1, 1);
    } else {
      k++;
    }
  }

  return b;
};

// Rational number
function Fraction(num, denom) {
  if (num instanceof Fraction)
    return num;

  if (typeof denom === 'object') {
    var g = num.gcd(denom);
    if (g.cmpn(1) !== 0) {
      num = num.div(g);
      denom = denom.div(g);
    }
  } else {
    denom = new bn(1);
  }

  // Ensure that denom is always positive
  if (denom.sign)
    num = num.neg();
  this.num = num;
  this.denom = denom;
}

Fraction.prototype.inspect = function inspect() {
  return '<Frac: ' +
      this.num.toString(10) + '/' + this.denom.toString(10) + '>';
};

Fraction.prototype.mul = function mul(f) {
  return new Fraction(this.num.mul(f.num), this.denom.mul(f.denom));
};

Fraction.prototype.sqr = function sqr() {
  return this.mul(this);
};

Fraction.prototype.sub = function sub(f) {
  // a/b - c/d = (a * d - c * b) / b * d
  var denom = this.denom.mul(f.denom);
  var num = this.num.mul(f.denom).isub(f.num.mul(this.denom));
  return new Fraction(num, denom);
};

Fraction.prototype.add = function add(f) {
  // a/b + c/d = (a * d + c * b) / b * d
  var denom = this.denom.mul(f.denom);
  var num = this.num.mul(f.denom).iadd(f.num.mul(this.denom));
  return new Fraction(num, denom);
};

Fraction.prototype.div = function div(f) {
  // (a/b) / (c/d) = a * d / b * c
  var num = this.num.mul(f.denom);
  var denom = this.denom.mul(f.num);
  return new Fraction(num, denom);
};

Fraction.prototype.cmpn = function cmpn(num) {
  var dm = this.num.divmod(this.denom);
  var r = dm.div.cmpn(num);
  if (r !== 0)
    return r;

  if (dm.mod.cmpn(0) !== 0)
    return this.num.sign ? 1 : -1;

  return 0;
};

Fraction.prototype.toBN = function toBN(round) {
  return round ? this.num.divRound(this.denom) : this.num.div(this.denom);
};

// Very primitive vector routines

function Vector(list) {
  this.length = list.length;
  this.list = new Array(this.length);
  for (var i = 0; i < this.length; i++)
    this.list[i] = new Fraction(list[i]);
}

Vector.prototype.inspect = function inspect() {
  return '<Vec: ' + this.list.map(function(item) {
    return item.inspect();
  }).join(', ') + '>';
};

Vector.prototype.dot = function dot(v) {
  var r = this.list[0].mul(v.list[0]);
  for (var i = 1; i < this.length; i++)
    r = r.add(this.list[i].mul(v.list[i]));
  return r;
};

Vector.prototype.sub = function sub(v) {
  var out = new Array(this.length);
  for (var i = 0; i < this.length; i++)
    out[i] = this.list[i].sub(v.list[i]);
  return new Vector(out);
};

Vector.prototype.mul = function mul(n) {
  var out = new Array(this.length);
  for (var i = 0; i < this.length; i++)
    out[i] = this.list[i].mul(n);
  return new Vector(out);
};
