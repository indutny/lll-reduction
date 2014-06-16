var assert = require('assert');

var lll = require('../');
var bn = require('bn.js');

describe('LLL', function() {
  describe('Gram-Schmidt process', function() {
    it('should make 2d-basis orthogonal', function() {
      var l = new lll();
      var basis = l.gramSchmidt([
        [ new bn(3), new bn(1) ],
        [ new bn(2), new bn(2) ]
      ]).basis;

      assert.equal(basis[0].dot(basis[1]).cmpn(0), 0);
    });

    it('should make 3d-basis orthogonal', function() {
      var l = new lll();
      var basis = l.gramSchmidt([
        [ new bn(3), new bn(2), new bn(1) ],
        [ new bn(1), new bn(3), new bn(2) ],
        [ new bn(2), new bn(1), new bn(3) ],
      ]).basis;

      assert.equal(basis[0].dot(basis[1]).cmpn(0), 0);
      assert.equal(basis[0].dot(basis[2]).cmpn(0), 0);
      assert.equal(basis[1].dot(basis[2]).cmpn(0), 0);
    });
  });

  describe('LLL itself', function() {
    it('should reduce 3d-basis', function() {
      var l = new lll();
      var basis = l.reduce([
        [ new bn(1), new bn(1), new bn(1) ],
        [ new bn(-1), new bn(0), new bn(2) ],
        [ new bn(3), new bn(5), new bn(6) ]
      ], new bn(3), new bn(4));

      var b = basis.map(function(b) {
        return b.inspect();
      }).join(', ');
      assert.equal(b, '<Vec: <Frac: 0/1>, <Frac: 1/1>, <Frac: 0/1>>, ' +
                      '<Vec: <Frac: 1/1>, <Frac: 0/1>, <Frac: 1/1>>, ' +
                      '<Vec: <Frac: -2/1>, <Frac: 0/1>, <Frac: 1/1>>');
    });

    it('should reduce 3d-basis', function() {
      var l = new lll();
      var basis = l.reduce([
        [ new bn(1), new bn(1), new bn(1) ],
        [ new bn(-1), new bn(0), new bn(2) ],
        [ new bn(3), new bn(5), new bn(6) ]
      ], new bn(3), new bn(4));

      var b = basis.map(function(b) {
        return b.inspect();
      }).join(', ');
      assert.equal(b, '<Vec: <Frac: 0/1>, <Frac: 1/1>, <Frac: 0/1>>, ' +
                      '<Vec: <Frac: 1/1>, <Frac: 0/1>, <Frac: 1/1>>, ' +
                      '<Vec: <Frac: -2/1>, <Frac: 0/1>, <Frac: 1/1>>');
    });
  });
});
