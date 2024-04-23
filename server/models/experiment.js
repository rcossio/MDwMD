const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const experimentSchema = new mongoose.Schema({
    diffusionCoefficient: Number,
    diffusionError: Number,
    diffusionUnit: String,
    manualProcessing: String,
    otherManualProcessing: String,
    method: String,
    otherMethod: String,
    proteinName: String,
    species: String,
    accessionNumber: String,
    referenceId: String,
    referenceIdType: String,
    referenceType: String,
    comment: String,
    timestamp: String,
    active: Boolean,
    userId: { type: mongoose.Schema.Types.ObjectId, ref: 'users' },
    supersededBy: { type: mongoose.Schema.Types.ObjectId, ref: 'experiments'},
});
  
// Create Mongoose Model
const ExperimentData = mongoose.model('experiments', experimentSchema);

module.exports = ExperimentData;
  