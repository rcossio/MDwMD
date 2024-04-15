const express = require('express');
const router = express.Router();
const verifyToken = require('../middleware/verifyToken');
const path = require('path');

//Public routes
router.get('/', (req, res) => {
  res.sendFile(path.resolve('public', 'browse.html'));
});

router.get('/login', (req, res) => {
  res.sendFile(path.resolve('public', 'login.html'));
});

//router.get('/register', (req, res) => {
//  res.sendFile(path.resolve('public', 'register.html'));
//});

//Protected routes
router.get('/new_experiment',verifyToken, (req, res) => {
  res.sendFile(path.resolve('public', 'new_experiment.html'));
});

router.get('/new_protein',verifyToken, (req, res) => {
  res.sendFile(path.resolve('public', 'new_protein.html'));
});

router.get('/manage/:id',verifyToken, (req, res) => {
  res.redirect('/manage.html?id='+req.params.id);
});

module.exports = router;
