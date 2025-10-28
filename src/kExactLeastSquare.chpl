module kExactLeastSquare 
{

class kls {
    var stencil_dom_ : domain(1) = {1..0};
    var stencilCellsCx_ : [stencil_dom_] real(64);
    var stencilCellsCy_ : [stencil_dom_] real(64);
    var x_ : [stencil_dom_] real(64);
    var y_ : [stencil_dom_] real(64);
    var r_ : [stencil_dom_] real(64);
    var c0_ : [stencil_dom_] real(64);
    var cx_ : [stencil_dom_] real(64);
    var cy_ : [stencil_dom_] real(64);

    var rho_fieldValues_ : [stencil_dom_] real(64);
    var u_fieldValues_ : [stencil_dom_] real(64);
    var v_fieldValues_ : [stencil_dom_] real(64);
    var E_fieldValues_ : [stencil_dom_] real(64);
    var p_fieldValues_ : [stencil_dom_] real(64);

    var centerCx_ : real(64);
    var centerCy_ : real(64);
    var A_ : [0..2, 0..2] real(64);
    var Ainv_ : [0..2, 0..2] real(64);


    var value_ : real(64);
    var gradX_ : real(64);
    var gradY_ : real(64);

    proc init(stencilCellsCx : [] real(64), stencilCellsCy : [] real(64), centerCx : real(64), centerCy : real(64)) {
        this.stencil_dom_ = stencilCellsCx.domain;
        this.stencilCellsCx_ = stencilCellsCx;
        this.stencilCellsCy_ = stencilCellsCy;
        this.centerCx_ = centerCx;
        this.centerCy_ = centerCy;
    }

    proc computeCoefficients() {
        for i in this.stencil_dom_ {
            const dx = this.stencilCellsCx_[i] - this.centerCx_;
            const dy = this.stencilCellsCy_[i] - this.centerCy_;
            const ri = 1 / sqrt(dx*dx + dy*dy);
            this.A_[0,0] += ri;
            this.A_[0,1] += ri * dx;
            this.A_[0,2] += ri * dy;
            this.A_[1,0] += ri * dx;
            this.A_[1,1] += ri * dx * dx;
            this.A_[1,2] += ri * dx * dy;
            this.A_[2,0] += ri * dy;
            this.A_[2,1] += ri * dx * dy;
            this.A_[2,2] += ri * dy * dy;

            this.x_[i] = dx;
            this.y_[i] = dy;
            this.r_[i] = ri;
        }

        // Solve for Ainv
        const detA = this.A_[0,0]*(this.A_[1,1]*this.A_[2,2] - this.A_[1,2]*this.A_[2,1]) -
                     this.A_[0,1]*(this.A_[1,0]*this.A_[2,2] - this.A_[1,2]*this.A_[2,0]) +
                     this.A_[0,2]*(this.A_[1,0]*this.A_[2,1] - this.A_[1,1]*this.A_[2,0]);
        this.Ainv_[0,0] =  (this.A_[1,1]*this.A_[2,2] - this.A_[1,2]*this.A_[2,1]) / detA;
        this.Ainv_[0,1] = -(this.A_[0,1]*this.A_[2,2] - this.A_[0,2]*this.A_[2,1]) / detA;
        this.Ainv_[0,2] =  (this.A_[0,1]*this.A_[1,2] - this.A_[0,2]*this.A_[1,1]) / detA;
        this.Ainv_[1,0] = -(this.A_[1,0]*this.A_[2,2] - this.A_[1,2]*this.A_[2,0]) / detA;
        this.Ainv_[1,1] =  (this.A_[0,0]*this.A_[2,2] - this.A_[0,2]*this.A_[2,0]) / detA;
        this.Ainv_[1,2] = -(this.A_[0,0]*this.A_[1,2] - this.A_[0,2]*this.A_[1,0]) / detA;
        this.Ainv_[2,0] =  (this.A_[1,0]*this.A_[2,1] - this.A_[1,1]*this.A_[2,0]) / detA;
        this.Ainv_[2,1] = -(this.A_[0,0]*this.A_[2,1] - this.A_[0,1]*this.A_[2,0]) / detA;
        this.Ainv_[2,2] =  (this.A_[0,0]*this.A_[1,1] - this.A_[0,1]*this.A_[1,0]) / detA;

        for i in this.stencil_dom_ {
            this.c0_[i] = this.r_[i] * (this.Ainv_[0,0] + this.Ainv_[0,1]*this.x_[i] + this.Ainv_[0,2]*this.y_[i]);
            this.cx_[i] = this.r_[i] * (this.Ainv_[1,0] + this.Ainv_[1,1]*this.x_[i] + this.Ainv_[1,2]*this.y_[i]);
            this.cy_[i] = this.r_[i] * (this.Ainv_[2,0] + this.Ainv_[2,1]*this.x_[i] + this.Ainv_[2,2]*this.y_[i]);
        }
    }

    proc interpolate(fieldValues : [] real(64)) {
        this.value_ = 0.0;
        this.gradX_ = 0.0;
        this.gradY_ = 0.0;
        for i in this.stencil_dom_ {
            this.value_ += this.c0_[i] * fieldValues[i];
            this.gradX_ += this.cx_[i] * fieldValues[i];
            this.gradY_ += this.cy_[i] * fieldValues[i];
        }

        return this.value_;
    }

}

}