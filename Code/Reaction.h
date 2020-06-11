//
// Created by Jacob Bradley on 2019-07-24.
//

#ifndef BIRTHSIMULATIONV1_REACTION_H
#define BIRTHSIMULATIONV1_REACTION_H

#include <iostream>
#include <random>
#include <forward_list>
using namespace std;

class Reaction {
private:
    int reactionType; // This will act as an index to track the type of reaction
    double completionTime; // This is the 'absolute system time' the reaction is scheduled to occur
public:
    // Default constructor
    Reaction(){
        reactionType = -1;
        completionTime = -1;
    }

    //TODO: add in paramiter checks - i.e. cant be negative
    Reaction(int tempReactionType, double tempCompletionTime){
        reactionType = tempReactionType;
        completionTime = tempCompletionTime;
    }

    int getReactionType(){
        return reactionType;
    }
    double getCompleationTime(){
        return  completionTime;
    }

    void setReactionType(int index){
        reactionType = index;
    }
    void setCompleationTime(double serviceTime, double systemTime, double tau){
        completionTime = serviceTime + systemTime + tau;
    }

    //TODO: How do destructors work again??
    //Destructor
    virtual ~Reaction() = default;

};

#endif //BIRTHSIMULATIONV1_REACTION_H
