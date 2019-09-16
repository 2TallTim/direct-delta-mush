#include <direct_delta_mush.h>

#include <maya/MArrayDataBuilder.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MItGeometry.h>
#include <maya/MMatrixArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MTimer.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#define EIGEN_NO_DEBUG

const MTypeId directDeltaMush::id(0x00131200);

void* directDeltaMush::creator()
{
    return new directDeltaMush();
}

MObject directDeltaMush::aTransSmooth;
// MObject directDeltaMush::aTransSmoothMap;
MObject directDeltaMush::aRotSmooth;
// MObject directDeltaMush::aRotSmoothMap;
MObject directDeltaMush::aAlpha;
MObject directDeltaMush::aRecalc;
MObject directDeltaMush::aSteps;

MStatus directDeltaMush::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnCompoundAttribute compAttr;

    MStatus status;

    aTransSmooth = nAttr.create("translationSmoothing", "ts", MFnNumericData::kDouble, 1.0, &status);
    nAttr.setMin(0.0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    addAttribute(aTransSmooth);

    // aTransSmoothMap = nAttr.create("translationSmoothingMap", "tsmap", MFnNumericData::kDouble, 1.0, &status);
    // nAttr.setMin(0.0);
    // nAttr.setMax(1.0);
    // nAttr.setArray(true);
    // nAttr.setUsesArrayDataBuilder(true);
    // addAttribute(aTransSmoothMap);

    aRotSmooth = nAttr.create("rotationSmoothing", "rs", MFnNumericData::kDouble, 1.0, &status);
    nAttr.setMin(0.0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    addAttribute(aRotSmooth);

    // aRotSmoothMap = nAttr.create("rotationSmoothingMap", "rsmap", MFnNumericData::kDouble, 1.0, &status);
    // nAttr.setMin(0.0);
    // nAttr.setMax(1.0);
    // nAttr.setArray(true);
    // nAttr.setUsesArrayDataBuilder(true);
    // addAttribute(aRotSmoothMap);

    aAlpha = nAttr.create("alpha", "alpha", MFnNumericData::kDouble);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    addAttribute(aAlpha);

    aSteps = nAttr.create("steps", "st", MFnNumericData::kInt, 1);
    nAttr.setMin(0);
    nAttr.setStorable(true);
    nAttr.setKeyable(true);
    addAttribute(aSteps);

    aRecalc = nAttr.create("recalc", "rc", MFnNumericData::kBoolean, false, &status);
    nAttr.setKeyable(true);
    nAttr.setStorable(true);
    addAttribute(aRecalc);

    attributeAffects(aTransSmooth, outputGeom);
    // attributeAffects(aTransSmoothMap, outputGeom);

    attributeAffects(aRotSmooth, outputGeom);
    // attributeAffects(aRotSmoothMap, outputGeom);

    attributeAffects(aAlpha, outputGeom);
    attributeAffects(aSteps, outputGeom);
    attributeAffects(aRecalc, outputGeom);

    // MGlobal::executeCommand("makePaintable -attrType multiDouble -sm deformer directDeltaMush translationSmoothingMap");
    // MGlobal::executeCommand("makePaintable -attrType multiDouble -sm deformer directDeltaMush rotationSmoothingMap");

    return MStatus::kSuccess;
}

void mmatrix_to_eigen(const MMatrix& M, Mat4& D)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            D.coeffRef(i, j) = M.matrix[i][j];
        }
    }
}

MMatrix eigen_to_mmatrix(const Mat4& E)
{
    double res[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res[i][j] = E.coeff(i, j);
        }
    }
    return MMatrix(res);
}

MStatus directDeltaMush::JumpToElement(MArrayDataHandle& hArray, unsigned int index)
{
    // Borrowed from Chad Vernon. Thanks Chad!
    MStatus status;
    status = hArray.jumpToElement(index);
    if (MFAIL(status)) {
        MArrayDataBuilder builder = hArray.builder(&status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        builder.addElement(index, &status);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = hArray.set(builder);
        CHECK_MSTATUS_AND_RETURN_IT(status);
        status = hArray.jumpToElement(index);
        CHECK_MSTATUS_AND_RETURN_IT(status);
    }
    return status;
}

MStatus
directDeltaMush::deform(MDataBlock& block,
    MItGeometry& iter,
    const MMatrix& /*m*/,
    unsigned int multiIndex)
{
    //
    // Method: deform
    //
    // Description:   Deforms the point with direct delta mush algorithm
    //
    // Arguments:
    //   block      : the datablock of the node
    //   iter       : an iterator for the geometry to be deformed
    //   m          : matrix to transform the point into world space
    //   multiIndex : the index of the geometry that we are deforming

    MStatus returnStatus;
    if (block.inputValue(aRecalc).asBool()) {
        precompute(block, iter);
    }

    if (omegas.empty()) {
        // Force precompute on load or create
        precompute(block, iter);
    }

    // get the influence transforms
    //
    MArrayDataHandle transformsHandle = block.inputArrayValue(matrix);
    int numTransforms = transformsHandle.elementCount();
    if (numTransforms == 0) {
        return MS::kSuccess;
    }

    MMatrixArray transforms;
    for (int i = 0; i < numTransforms; ++i) {
        transforms.append(MFnMatrixData(transformsHandle.inputValue().data()).matrix());
        transformsHandle.next();
    }

    MArrayDataHandle bindHandle = block.inputArrayValue(bindPreMatrix);
    if (bindHandle.elementCount() > 0) {
        for (auto& transform : transforms) {
            MMatrix pb = MFnMatrixData(bindHandle.inputValue().data()).matrix();
            transform = pb * transform;
            bindHandle.next();
        }
    }

    MArrayDataHandle weightListHandle = block.inputArrayValue(weightList);
    if (weightListHandle.elementCount() == 0) {
        // no weights - nothing to do
        return MS::kSuccess;
    }

    // Iterate through each point in the geometry.
    //
    MPointArray positions;
    iter.allPositions(positions);

    for (int idx = 0; idx < positions.length(); idx++) {
        MPoint& pt = positions[idx];
        Mat4 qmat;
        qmat.setZero();
        for (auto e : omegas.at(idx)) {
            Mat4 Mi;
            mmatrix_to_eigen(transforms[e.first], Mi);
            Mat4 tmp = (getOmega(idx, e.first) * Mi);
            qmat += tmp;
        }

        qmat.transposeInPlace();

        Mat3 Qi = qmat.block(0, 0, 3, 3);
        Vec3 qi = qmat.block(0, 3, 3, 1);
        Vec3 pi = qmat.block(3, 0, 1, 3).transpose();

        // Rotation
        Mat3 M = Qi - (qi * pi.transpose());
        Eigen::JacobiSVD<Mat3> svd;
        svd.compute(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Mat3 U = svd.matrixU();
        Mat3 V = svd.matrixV().transpose();
        Mat3 R = (U * V);

        //Translation
        Vec3 ti = qi - (R * pi);

        Mat4 gamma;
        gamma << R, ti, 0, 0, 0, 1;

        Vec4 pt_h(pt.x, pt.y, pt.z, 1);

        Vec4 final_pt = gamma * pt_h;

        pt = MPoint(final_pt[0], final_pt[1], final_pt[2]);
    }
    iter.setAllPositions(positions);
    return returnStatus;
}

// standard initialization procedures
//

MStatus initializePlugin(MObject obj)
{
    MStatus result;

    MFnPlugin plugin(obj, "Tim Henning", "0.0.1", "Any");
    result = plugin.registerNode(
        "directDeltaMush",
        directDeltaMush::id,
        &directDeltaMush::creator,
        &directDeltaMush::initialize,
        MPxNode::kSkinCluster);

    return result;
}

Vec4 hmg(const Vec3& v)
{
    auto res = Vec4(v[0], v[1], v[2], 1);
    return res;
}

Mat4 directDeltaMush::getOmega(int i, int j)
{
    //Get omega, or if it doesn't exist, return the zero matrix.
    auto l = omegas.at(i);
    for (auto e : l) {
        if (e.first == j)
            return e.second;
    }
    Mat4 z;
    z.setZero();
    return z;
}

MStatus directDeltaMush::precompute(MDataBlock& block, MItGeometry& iter)
{
    // Time the precomputation stage
    MTimer timer;
    timer.beginTimer();
    MGlobal::displayInfo("Starting Precomputation");

    // Get the geometry
    MArrayDataHandle hInput = block.inputArrayValue(input);
    JumpToElement(hInput, 0);
    MDataHandle hInputGeom = hInput.inputValue().child(inputGeom);
    MObject inmesh = hInputGeom.asMesh();
    MFnMesh mesh(inmesh);
    // Verts
    Eigen::MatrixX3d V; // n by 3 matrix of undeformed vertex positions
    MPointArray verts;
    mesh.getPoints(verts);
    V.resize(verts.length(), Eigen::NoChange);
    for (unsigned int i = 0; i < verts.length(); i++) {
        double vx = verts[i][0];
        double vy = verts[i][1];
        double vz = verts[i][2];
        V.row(i) = Eigen::RowVector3d(vx, vy, vz);
    }

    // Faces
    MIntArray tricount, triidx;
    Eigen::MatrixX3i F; // n by 3 matrix of undeformed triangle indices
    mesh.getTriangles(tricount, triidx);
    F.resize(triidx.length() / 3, Eigen::NoChange);
    for (unsigned int i = 0; i < triidx.length() / 3; i++) {
        F.row(i) = Eigen::RowVector3i(triidx[(3 * i)], triidx[(3 * i) + 1], triidx[(3 * i) + 2]);
    }

    std::ptrdiff_t n = V.rows(); // Number of vertices

    // Transform matrices

    MArrayDataHandle transformsHandle = block.inputArrayValue(matrix);
    int numTransforms = transformsHandle.elementCount();
    if (numTransforms == 0) {
        return MS::kSuccess;
    }

    MMatrixArray transforms;
    for (int i = 0; i < numTransforms; ++i) {
        transforms.append(MFnMatrixData(transformsHandle.inputValue().data()).matrix());
        transformsHandle.next();
    }

    // Weights & other maps
    Sparse W(n, numTransforms);
    W.reserve(4); // Assuming approximately 4 weights per vertex

    //Eigen::VectorXd rotSmoothMap(n);
    //Eigen::VectorXd transSmoothMap(n);

    MArrayDataHandle hWeightList = block.inputArrayValue(weightList);
    //MArrayDataHandle hRotSmoothMap = block.inputArrayValue(aRotSmoothMap);
    //MArrayDataHandle hTransSmoothMap = block.inputArrayValue(aTransSmoothMap);

    long ii = 0;
    for (iter.reset(); !iter.isDone(); iter.next(), ii++) {
        MArrayDataHandle weightsHandle = hWeightList.inputValue().child(weights);
        for (long widx = 0; widx < numTransforms; widx++) {
            JumpToElement(weightsHandle, widx);
            double w = weightsHandle.inputValue().asDouble();
            if (w != 0 && !isnan(w)) {
                W.insert(ii, widx) = w;
            }
        }
        //JumpToElement(hRotSmoothMap, iter.index());
        //rotSmoothMap[ii] = hRotSmoothMap.inputValue().asFloat();
        //JumpToElement(hTransSmoothMap, iter.index());
        //transSmoothMap[ii] = hTransSmoothMap.inputValue().asFloat();

        hWeightList.next();
    }

    // Get precomputation parameters
    double rotSmooth = block.inputValue(aRotSmooth).asDouble();
    double transSmooth = block.inputValue(aTransSmooth).asDouble();

    double dmBlend = block.inputValue(aAlpha).asDouble();

    int steps = block.inputValue(aSteps).asInt();

    //Laplacian matrix
    Sparse lapl;
    igl::cotmatrix(V, F, lapl); // Compute standard laplacian matrix
    MatX lapl_diag_inv = lapl.diagonal().asDiagonal().inverse(); //Normalize
    Sparse L = (lapl * lapl_diag_inv).sparseView().eval(); // Normalized!

    // Vars needed for solver. Solver is used to calculate sparse inverses.
    Sparse I(n, n);
    I.setIdentity();

    Eigen::SparseLU<Sparse> solver_b;
    Eigen::SparseLU<Sparse> solver_c;

    // Implicitly solve.
    // This is a slight deviation from the original paper, which has parameters for
    // steps and the separate translation and rotation smoothing. For an artist, it's
    // easier to think in terms of the total desired amount of smoothing and the
    // number of steps as a precision parameter, rather than having to tune them in
    // tandem for each change.
    Sparse b(I + (transSmooth / (double)steps) * L /*transSmoothMap.asDiagonal()*/);
    Sparse c(I + (rotSmooth / (double)steps) * L /*rotSmoothMap.asDiagonal()*/);

    Sparse B(I);
    Sparse B_next;
    solver_b.compute(b.transpose());
    for (int i = 0; i < steps; i++) {
        B.makeCompressed();
        // This is the slow part
        B_next = solver_b.solve(B);
        B = B_next;
    }

    Sparse C(I);
    Sparse C_next;
    solver_c.compute(c.transpose());
    for (int i = 0; i < steps; i++) {
        C.makeCompressed();
        C_next = solver_c.solve(C);
        C = C_next;
    }

    auto psi = [&B, &V, &W, n, &transforms](const long i, const long j) -> Mat4 {
        Mat4 sum = Mat4::Zero();
        // Sum psi, a sum over all verts for a given weight.

        for (long k = 0; k < n; k++) {
            double w = W.coeff(k, j);
            double b = B.coeff(k, i);
            if (w != 0 && b != 0) {
                Vec3 r = V.row(k);
                Vec4 rw(V.coeff(k, 0), V.coeff(k, 1), V.coeff(k, 2), 1);
                Mat4 hh = (rw * rw.transpose());
                Mat4 h = b * w * hh;
                sum += h;
                if (h != h.transpose()) {
                    h = h.transpose();
                }
            }
        }
        return sum;
    };

    // Helper lambdas for the precomputation process.

    auto p_i = [&psi, numTransforms](const int ii) -> Vec3 {
        Mat4 sum = Mat4::Zero();
        for (int j = 0; j < numTransforms; j++) {
            sum += psi(ii, j);
        }
        if (sum != sum.transpose()) {
            throw sum;
        }
        return sum.block(0, 3, 3, 1).eval();
    };

    auto w_prime = [&C, &W, n](const int i, const int j) -> double {
        double sum = 0;
        for (int k = 0; k < n; k++) {
            double w = W.coeff(k, j);
            double c = C.coeff(k, i);
            sum += w * c;
        }
        return sum;
    };

    auto omega = [dmBlend, &psi, &p_i, &w_prime, &W](const int i, const int j) -> Mat4 {
        Vec3 p_ii = p_i(i);
        Mat4 pimat;
        pimat << (p_ii * p_ii.transpose()), p_ii, p_ii.transpose(), 1;
        Mat4 psi_ij = psi(i, j);
        return ((1.0 - dmBlend) * (psi_ij) + (dmBlend * w_prime(i, j)) * pimat).eval();
    };

    omegas.clear();
    // Actually precompute omegas.
    for (int ii = 0; ii < n; ii++) {
        omegas.push_back(std::list<std::pair<int, Mat4>>());
        auto& e = omegas.at(ii);
        e.clear();
        for (int jj = 0; jj < numTransforms; jj++) {
            // This could be optimized more by not storing zero matrices
            //if(W.coeff(ii,jj) != 0){
            auto o = omega(ii, jj);
            e.push_back(std::pair<int, Mat4>(jj, o));
            //}
        }
    }

    timer.endTimer();

    std::stringstream ss;
    ss << "Precomputation complete (" << timer.elapsedTime() << ")";
    MGlobal::displayInfo(ss.str().c_str());
    return MS::kSuccess;
}

MStatus uninitializePlugin(MObject obj)
{
    MStatus result;

    MFnPlugin plugin(obj);
    result = plugin.deregisterNode(directDeltaMush::id);

    return result;
}
