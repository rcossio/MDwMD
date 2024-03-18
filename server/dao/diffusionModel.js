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
    userName: String,
    active: Boolean,
    supersededBy: { type: mongoose.Schema.Types.ObjectId },
    pdbStructures: [String],
    discardedPdbStructures: [{ pdbId: String, reason: String, _id: false }]
});
  
// Create Mongoose Model
const DiffusionData = mongoose.model('data_from_experiments', diffusionSchema);

module.exports = DiffusionData;
  