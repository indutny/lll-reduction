var assert = require('assert');

var lll = require('../');
var F = require('fract.js');

describe('LLL', function() {
  describe('Gram-Schmidt process', function() {
    it('should make 2d-basis orthogonal', function() {
      var l = new lll();
      var basis = l.gramSchmidt([
        [ new F(3), new F(1) ],
        [ new F(2), new F(2) ]
      ]).basis;

      assert.equal(basis[0].dot(basis[1]).cmpn(0), 0);
    });

    it('should make 3d-basis orthogonal', function() {
      var l = new lll();
      var basis = l.gramSchmidt([
        [ new F(3), new F(2), new F(1) ],
        [ new F(1), new F(3), new F(2) ],
        [ new F(2), new F(1), new F(3) ],
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
        [ new F(1), new F(1), new F(1) ],
        [ new F(-1), new F(0), new F(2) ],
        [ new F(3), new F(5), new F(6) ]
      ], new F(3, 4));

      var b = basis.map(function(b) {
        return b.inspect();
      }).join(', ');
      assert.equal(b, '<Vec: <Frac: 0/1>, <Frac: 1/1>, <Frac: 0/1>>, ' +
                      '<Vec: <Frac: 1/1>, <Frac: 0/1>, <Frac: 1/1>>, ' +
                      '<Vec: <Frac: -2/1>, <Frac: 0/1>, <Frac: 1/1>>');
    });
  });
});
