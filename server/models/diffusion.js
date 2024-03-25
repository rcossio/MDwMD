const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const diffusionSchema = new mongoose.Schema({
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
    pdbStructures: [String],
    discardedPdbStructures: [{ pdbId: String, reason: String, _id: false }]
});
  
// Create Mongoose Model
const DiffusionData = mongoose.model('experiments', diffusionSchema);

module.exports = DiffusionData;
  