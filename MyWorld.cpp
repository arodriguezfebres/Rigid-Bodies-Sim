#include "MyWorld.h"
#include "RigidBody.h"
#include "CollisionInterface.h"
#include <iostream>

using namespace Eigen;
using namespace std;

MyWorld::MyWorld() {
    mFrame = 0;
    mTimeStep = 0.001;
    mGravity = Vector3d(0.0, -9.8, 0.0);
    mForce.setZero();
    // Create a collision detector
    mCollisionDetector = new CollisionInterface();
    
    // Create and intialize two default rigid bodies
    RigidBody *rb1 = new RigidBody(dart::dynamics::Shape::BOX, Vector3d(0.05, 0.05, 0.05));
    mCollisionDetector->addRigidBody(rb1, "box"); // Put rb1 in collision detector
    rb1->mPosition[0] = -0.3;
    rb1->mPosition[1] = -0.5;
    
    rb1->mAngMomentum = Vector3d(0.0, 0.01, 0.0);
    mRigidBodies.push_back(rb1);
    
    RigidBody *rb2 = new RigidBody(dart::dynamics::Shape::ELLIPSOID, Vector3d(0.06, 0.06, 0.06));
    mCollisionDetector->addRigidBody(rb2, "ellipse"); // Put rb2 in collision detector
    rb2->mPosition[0] = 0.3;
    rb2->mPosition[1] = -0.5;
    rb2->mAngMomentum = Vector3d(0.01, 0.0, 0.0);
    rb2->mColor = Vector4d(0.2, 0.8, 0.2, 1.0); // Blue
    mRigidBodies.push_back(rb2);
}

void MyWorld::initializePinata() {
    // Add pinata to the collison detector
    mCollisionDetector->addSkeleton(mPinataWorld->getSkeleton(0));
    
    // Add some damping in the Pinata joints
    int nJoints = mPinataWorld->getSkeleton(0)->getNumBodyNodes();
    for (int i = 0; i < nJoints; i++) {
        int nDofs = mPinataWorld->getSkeleton(0)->getJoint(i)->getNumDofs();
        for (int j = 0; j < nDofs; j++)
        mPinataWorld->getSkeleton(0)->getJoint(i)->setDampingCoefficient(j, 1.0);
    }
    
    // Weld two seems to make a box
    dart::dynamics::BodyNode* top = mPinataWorld->getSkeleton(0)->getBodyNode("top");
    dart::dynamics::BodyNode* front = mPinataWorld->getSkeleton(0)->getBodyNode("front");
    dart::dynamics::BodyNode* back = mPinataWorld->getSkeleton(0)->getBodyNode("back");
    dart::constraint::WeldJointConstraint *joint1 = new dart::constraint::WeldJointConstraint(top, front);
    dart::constraint::WeldJointConstraint *joint2 = new dart::constraint::WeldJointConstraint(top, back);
    mPinataWorld->getConstraintSolver()->addConstraint(joint1);
    mPinataWorld->getConstraintSolver()->addConstraint(joint2);
}

MyWorld::~MyWorld() {
    for (int i = 0; i < mRigidBodies.size(); i++)
    delete mRigidBodies[i];
    mRigidBodies.clear();
    if (mCollisionDetector)
    delete mCollisionDetector;
}

void MyWorld::simulate() {
    mFrame++;
    
	MatrixXd iBody = mRigidBodies[0]->mShape->computeInertia(mRigidBodies[0]->mMass);
	MatrixXd iTensor;
	MatrixXd angularVelocity;

    // TODO: The skeleton code has provided the integration of position and linear momentum,
    // your first job is to fill in the integration of orientation and angular momentum.
    for (int i = 0; i < mRigidBodies.size(); i++) {
        // derivative of position and linear momentum
        Eigen::Vector3d dPos = mRigidBodies[i]->mLinMomentum / mRigidBodies[i]->mMass;
        Eigen::Vector3d dLinMom = mRigidBodies[i]->mMass * mGravity + mRigidBodies[i]->mAccumulatedForce;
        
        // update position and linear momentum
        mRigidBodies[i]->mPosition += dPos * mTimeStep;
        mRigidBodies[i]->mLinMomentum += mTimeStep * dLinMom;

		// Uhhhh? Solve for this stuff.
		MatrixXd tempRotMatrix = mRigidBodies[i]->mQuatOrient.toRotationMatrix();
		iTensor = tempRotMatrix * iBody * tempRotMatrix.transpose();
		angularVelocity = iTensor.inverse() * mRigidBodies[i]->mAngMomentum; //w?
		mRigidBodies[i]->mQuatOrient.normalize();


		//QDOT - START:
		Quaterniond aVel = Quaterniond(0.0, angularVelocity(0), angularVelocity(1), angularVelocity(2));
		
		aVel.w() = .5 * aVel.w();
		aVel.vec() = .5 * aVel.vec();
		Quaterniond tempHolder;
		Quaterniond crb = mRigidBodies[i]->mQuatOrient;

		tempHolder.w() = (aVel.w() * crb.w()) - (aVel.vec()(0) * crb.vec()(0) + aVel.vec()(1) * crb.vec()(1) + aVel.vec()(2) * crb.vec()(2));
		Vector3d crossCut = Vector3d(aVel.vec()(0), aVel.vec()(1), aVel.vec()(2)).cross(Vector3d(crb.vec()(0), crb.vec()(1), crb.vec()(2)));
		crossCut = crossCut + (aVel.w() * crb.vec()) + (crb.w() * aVel.vec());
		tempHolder.vec()(0) = crossCut(0);
		tempHolder.vec()(1) = crossCut(1);
		tempHolder.vec()(2) = crossCut(2);

		tempHolder.w() = tempHolder.w() * mTimeStep;
		tempHolder.vec() = tempHolder.vec() * mTimeStep;

		//QDOT - End

		mRigidBodies[i]->mQuatOrient.w() = mRigidBodies[i]->mQuatOrient.w() + tempHolder.w();
		mRigidBodies[i]->mQuatOrient.vec() = mRigidBodies[i]->mQuatOrient.vec() + tempHolder.vec();
		mRigidBodies[i]->mAngMomentum = mRigidBodies[i]->mAngMomentum + mRigidBodies[i]->mAccumulatedTorque * mTimeStep;

		//update orientation
		mRigidBodies[i]->mOrientation = mRigidBodies[i]->mQuatOrient.toRotationMatrix();
	}
    
    // Reset accumulated force and torque to be zero after a complete integration
    for (int i = 0; i < mRigidBodies.size(); i++) {
        mRigidBodies[i]->mAccumulatedForce.setZero();
        mRigidBodies[i]->mAccumulatedTorque.setZero();
    }
    
    // Apply external force to the pinata
    mPinataWorld->getSkeleton(0)->getBodyNode("bottom")->addExtForce(mForce);
    mForce.setZero();
    
    // Simulate Pinata using DART
    mPinataWorld->step();
    
    // Run collision detector
    mCollisionDetector->checkCollision();
    
    // TODO: implement a collision handler
    collisionHandling();
    
    // Break the pinata if it has enough momentum
    if (mPinataWorld->getSkeleton(0)->getCOMLinearVelocity().norm() > 0.6)
    mPinataWorld->getConstraintSolver()->removeAllConstraints();
}

// TODO: fill in the collision handling function
void MyWorld::collisionHandling() {
    // restitution coefficient
    double epsilon = 0.8;
    
    // TODO: handle the collision events
	int contacts = mCollisionDetector->getNumContacts();

	for (int index = 0; index < contacts; index++) {
		RigidBody* myA = mCollisionDetector->getContact(index).rb1;
		RigidBody* myB = mCollisionDetector->getContact(index).rb2;

		Vector3d myN = mCollisionDetector->getContact(index).normal;
		Vector3d pVel = mCollisionDetector->getContact(index).pinataVelocity;

		Vector3d aPDotMinus;
		Vector3d bPDotMinus;

		if (myA != nullptr && myB != nullptr) {
			MatrixXd aIBody = myA->mShape->computeInertia(myA->mMass);
			MatrixXd A_Inot = (myA->mQuatOrient.toRotationMatrix() * aIBody * myA->mQuatOrient.toRotationMatrix().transpose()).inverse();
			Vector3d A_w = A_Inot * myA->mAngMomentum;
			Vector3d A_r = mCollisionDetector->getContact(index).point - myA->mPosition;
			aPDotMinus = myA->mLinMomentum + (A_w.cross(A_r));
			
			MatrixXd bIBody = myB->mShape->computeInertia(myB->mMass);
			MatrixXd B_Inot = (myB->mQuatOrient.toRotationMatrix() * bIBody * myB->mQuatOrient.toRotationMatrix().transpose()).inverse();
			Vector3d B_w = B_Inot * myB->mAngMomentum;
			Vector3d B_r = mCollisionDetector->getContact(index).point - myB->mPosition;
			bPDotMinus = myB->mLinMomentum + (B_w.cross(B_r));

			float doWeAccept = myN.dot(aPDotMinus - bPDotMinus);

			if (doWeAccept < 0.0) {
				float myJ = ((-1 * (1 + epsilon)) * doWeAccept) / ((1.0 / myA->mMass) + (1.0 / myB->mMass) + (myN.dot((A_Inot * A_r.cross(myN).cross(A_r)))) + (myN.dot((B_Inot * B_r.cross(myN).cross(B_r)))));
				myA->mLinMomentum = myA->mLinMomentum + (myJ * myN);
				myA->mAngMomentum = myA->mAngMomentum + (A_r.cross(myJ * myN));

				myB->mLinMomentum = myB->mLinMomentum + (-myJ * myN);
				myB->mAngMomentum = myB->mAngMomentum + (B_r.cross(-myJ * myN));
			}

		} else if (myA != nullptr && myB == nullptr) {
			MatrixXd aIBody = myA->mShape->computeInertia(myA->mMass);
			MatrixXd A_Inot = (myA->mQuatOrient.toRotationMatrix() * aIBody * myA->mQuatOrient.toRotationMatrix().transpose()).inverse();
			Vector3d A_w = A_Inot * myA->mAngMomentum;
			Vector3d A_r = mCollisionDetector->getContact(index).point - myA->mPosition;
			aPDotMinus = myA->mLinMomentum + (A_w.cross(A_r));
			bPDotMinus = mCollisionDetector->getContact(index).pinataVelocity;

			float doWeAccept = myN.dot(aPDotMinus - bPDotMinus);

			if (doWeAccept < 0.0) {
				float myJ = ((-1 * (1 + epsilon)) * doWeAccept) / ((1.0 / myA->mMass) + (myN.dot((A_Inot * A_r.cross(myN).cross(A_r)))) );
				myA->mLinMomentum = myA->mLinMomentum + (myJ * myN);
				myA->mAngMomentum = myA->mAngMomentum + (A_r.cross(myJ * myN));
			}

		} else if (myA == nullptr && myB != nullptr) {
			MatrixXd bIBody = myB->mShape->computeInertia(myB->mMass);
			MatrixXd B_Inot = (myB->mQuatOrient.toRotationMatrix() * bIBody * myB->mQuatOrient.toRotationMatrix().transpose()).inverse();
			Vector3d B_w = B_Inot * myB->mAngMomentum;
			Vector3d B_r = mCollisionDetector->getContact(index).point - myB->mPosition;
			bPDotMinus = myB->mLinMomentum + (B_w.cross(B_r));
			aPDotMinus = mCollisionDetector->getContact(index).pinataVelocity;

			float doWeAccept = myN.dot(aPDotMinus - bPDotMinus);

			if (doWeAccept < 0.0) {
				float myJ = ((-1 * (1 + epsilon)) * doWeAccept) / ((1.0 / myB->mMass) + (myN.dot((B_Inot * B_r.cross(myN).cross(B_r)))));
				myB->mLinMomentum = myB->mLinMomentum + (-myJ * myN);
				myB->mAngMomentum = myB->mAngMomentum + (B_r.cross(-myJ * myN));
			}
		}


	}
}









